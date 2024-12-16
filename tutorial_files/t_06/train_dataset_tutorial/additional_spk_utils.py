import subprocess as sbp
import os
from ase import Atoms
import schnetpack as spk
from shutil import rmtree
from mopac import System_Properties
from glob import glob

class Run_Training():
    def __init__(self,
                 remove_prev_model: bool=True, 
                 db_path: str=None,
                 model_directory: str=None,
                 property: str='energy',
                 derivative: str='forces',
                 rho: tuple=(0.01, 0.99),
                 split: tuple=(None, None),
                 batch_size: int=8,
                 device: str='cuda',
                 n_epochs: int=5000,
                 lr: float=0.0001,
                 lr_patience: int=10,
                 lr_min: float=1e-06,
                 cutoff: float=100.0,
                 num_gaussians: int=50,
                 features: int=128,
                 interactions: int=3) -> None:
        self.remove_prev_model=remove_prev_model
        self.db_path=db_path
        self.model_dir=model_directory
        self.property=property
        self.derivative=derivative
        self.rho=rho
        self.split=split
        self.batch_size=batch_size
        self.device=device
        self.n_epochs=n_epochs
        self.lr=lr
        self.lr_patience=lr_patience
        self.lr_min=lr_min
        self.cutoff=cutoff
        self.num_gaussians=num_gaussians
        self.features=features
        self.interactions=interactions
        return None
    
    def create_training_script(self):
        with open('train.sh', mode='w', encoding='utf8') as w_file:
            w_file.write('#!/bin/bash\n')
            if self.device == 'cuda':
                w_file.write(f'python spk_run.py train schnet custom {self.db_path} {self.model_dir} --property {self.property} --derivative {self.derivative} --rho property={self.rho[0]} derivative={self.rho[1]} --split {self.split[0]} {self.split[1]} --batch_size {self.batch_size} --{self.device} --n_epochs {self.n_epochs} --lr {self.lr} --lr_patience {self.lr_patience} --lr_min {self.lr_min} --cutoff {self.cutoff} --num_gaussians {self.num_gaussians} --features {self.features} --interactions {self.interactions}\n')
                w_file.write(f'python spk_run.py eval {self.db_path} {self.model_dir} --split test --{self.device} --batch_size 1\n')
            else:
                w_file.write(f'python spk_run.py train schnet custom {self.db_path} {self.model_dir} --property {self.property} --derivative {self.derivative} --rho property={self.rho[0]} derivative={self.rho[1]} --split {self.split[0]} {self.split[1]} --batch_size {self.batch_size} --n_epochs {self.n_epochs} --lr {self.lr} --lr_patience {self.lr_patience} --lr_min {self.lr_min} --cutoff {self.cutoff} --num_gaussians {self.num_gaussians} --features {self.features} --interactions {self.interactions}\n')
                w_file.write(f'python spk_run.py eval {self.db_path} {self.model_dir} --split test --batch_size 1\n')
    
    def run(self):
        sbp.run(['chmod', 'u+x', 'train.sh'])
        model_path = os.path.join('.', self.model_dir)
        if self.remove_prev_model and os.path.isdir(model_path):
            rmtree(model_path)
        eval_file = os.path.join(model_path, 'evaluation.txt')
        if os.path.isfile(eval_file):
            os.remove(eval_file)
        res = sbp.run(['./train.sh'], shell=True, stderr=sbp.PIPE, stdout=sbp.PIPE, text=True)
        print(f'process finished with the return code: {res.returncode}')
        if res.returncode != 0:
            print(f'The process finished with the following error: {res.stderr}')

    def run_training(self):
        self.create_training_script()
        self.run()

class Build_AseDb():
    def __init__(self, load_existing_database: bool = False, db_name: str=None, db_properties: list=None, metadata: dict=None) -> None:
        self.load_existing_db = load_existing_database
        if db_name.count('.db') == 0:
            self.db_path = os.path.join('.', f'{db_name}.db')
        else:
            self.db_path = os.path.join('.', db_name)
        if os.path.isfile(self.db_path) and self.load_existing_db == False:
            print(f'Removing preivous database: {self.db_path}')
            os.remove(self.db_path)
        self.db_properties = db_properties
        self.metadata = metadata
        return None
    
    def create_db(self):
        if self.load_existing_db:
            try:
                return spk.data.AtomsData(self.db_path)
            except spk.data.AtomsDataError: print('The specified database does not exists!\nPlease specify correct name of the database.')
        else:
            db = spk.data.AtomsData(self.db_path, available_properties=self.db_properties)
            db.set_metadata(metadata=self.metadata)
            return db

class BuildAsedbError(Exception):
    pass

class Build_AseDb_From_Mopac_Aux(Build_AseDb):
    def __init__(self, load_existing_database: bool = False,
                 complex_path: str=None,
                 buffer_path: str=None,
                 db_name: str = None, 
                 db_properties: tuple = None, 
                 metadata: dict = None, 
                 inner_region_size: int=None, 
                 reference_energies: tuple=None) -> None:
        super().__init__(load_existing_database, db_name, db_properties, metadata)
        self.complex_path = complex_path
        self.complex_files = sorted(glob(os.path.join(self.complex_path, '*.aux')))
        self.buffer_path = buffer_path
        self.buffer_files = sorted(glob(os.path.join(self.buffer_path, '*.aux')))
        self.inner_region_size = inner_region_size
        self.reference_energies = reference_energies

    def check(self, files, path):
        if len(files) == 0:
            raise BuildAsedbError(f'No aux files found on the path: {path}\nMake sure that the directory MOPAC_results is in this directory!')

    def get_mopac_properties(self, complex_path, buffer_path)->dict:
        props = System_Properties(complex_path, buffer_path, self.inner_region_size, self.reference_energies)
        sys_props = {}
        sys_props['complex_energy'], sys_props['buffer_energy'] = props.get_energy()
        sys_props['num_of_waters_complex'], sys_props['num_of_waters_buffer'] = props.get_num_of_h2o()
        sys_props['complex_atoms'], sys_props['buffer_atoms'] = props.get_atoms()
        sys_props['complex_positions'], sys_props['buffer_positions'] = props.get_positions()
        sys_props['complex_forces'], sys_props['buffer_forces'] = props.get_forces()
        sys_props['energy'] = props.get_burnn_energy()
        sys_props['forces'] = props.get_burnn_forces()
        sys_props['spin'] = props.get_spin()
        return sys_props

    def build_db(self):
        self.check(self.complex_files, self.complex_path)
        self.check(self.buffer_files, self.buffer_path)
        db = self.create_db()
        for complex_path, buffer_path in zip(self.complex_files, self.buffer_files):
            system_props = self.get_mopac_properties(complex_path=complex_path, buffer_path=buffer_path)
            db_props = {}
            for item in self.db_properties:
                try:
                    db_props[item] = system_props[item]
                except KeyError: raise BuildAsedbError(f'Property {item} is not possible to get from MOPAC now.\nallowed properties are:\ncomplex_energy\nbuffer_energy\nnum_of_waters_complex\nnum_of_waters_buffer\ncomplex_atoms\nbuffer_atoms\ncomplex_positions\nbuffer_positions\ncomplex_forces\nbuffer_forces\nenergy\nforces\nspin')
            system = Atoms(symbols=system_props['complex_atoms'], positions=system_props['complex_positions'])
            db.add_system(system, properties=db_props)