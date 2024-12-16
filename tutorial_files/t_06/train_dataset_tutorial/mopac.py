import os, re
import logging as log
import numpy as np
import multiprocessing as mp
import subprocess as sbp
import logging.config as logcnf
import atomic_properties
from gromos import ReShake
from glob import glob
from shutil import rmtree

# os.system('rm *.log')

#logcnf.fileConfig('/home/user/Python/Python_modules/logging.info')
#logger = log.getLogger(__name__)

KCAL_TO_KJ = 4.184

class MopacCalculationError(Exception):
    pass

class WrongAtomTypeError(Exception):
    pass

class Mopac_Calculation():
    '''
    Basic class for Mopac calculations.
    The class creates file system and common
    class variables.

    Input: coord_file: array like object of configuration files
    Methods: find_atom_type: returns atomic number (required by Mopac)
             get_positions_in_A: convers GROMOS coordinates (nm) to Angstroms
             create_run_files: creates scripts for running individual Mopac calculations
             run_calculation: parallelization and execution of scripts created
                              by create_run_files method using method run
    
    Only Mopac input files, scripts for running the program and Mopac 
    AUX output files are stored. Other Mopac output files are removed.
    '''
    def __init__(self, coord_files: list, inner_region_size: int = None, 
                 topo_path: str = None, gromos_bin: str = None, imd_example: str = None,
                 mk_example: str = None, mk_lib: str = None):
        self.coord_files = coord_files
        self.mopac_dir = os.path.join(os.getcwd(), 'MOPAC_results')
        self.topo_path = topo_path
        self.bin = gromos_bin
        self.imd_example = imd_example
        self.mk_example = mk_example
        self.lib = mk_lib
        self.step = None

        # input from GROMOS
        if self.coord_files[0].endswith('.cnf'):
            try:
                os.mkdir(self.mopac_dir) # creation of the directory for the mopac files
                os.mkdir(os.path.join(self.mopac_dir, 'buffer_pls_inner'))
                os.mkdir(os.path.join(self.mopac_dir, 'buffer'))
            except FileExistsError: pass
            self.mopac_path = os.path.join(self.mopac_dir, 'buffer_pls_inner')
            #logger.info(f'Directory {self.mopac_path} was successfully created.')
            #logger.info(f'MOPAC path was set to {self.mopac_path}')

        # input from AUX file
        elif self.coord_files[0].endswith('.aux'):
            try:
                os.mkdir(os.path.join(self.mopac_dir, 'minimization')) # creation of the directory for the mopac files
                os.mkdir(os.path.join(self.mopac_dir, 'minimization', 'buffer_pls_inner'))
                os.mkdir(os.path.join(self.mopac_dir, 'minimization', 'buffer'))
            except FileExistsError: pass
            self.mopac_path = os.path.join(self.mopac_dir, 'minimization', 'buffer_pls_inner')
            #logger.info(f'Directory {self.mopac_path} was successfully created.')
        try:
            os.mkdir(f'{self.mopac_path}/aux_out')
            if self.mopac_path != os.path.join(self.mopac_dir, 'minimization', 'buffer_pls_inner'):
                os.mkdir(f'{self.mopac_path}/run_files')
            #logger.info(f'Directories {self.mopac_path}/run_files and {self.mopac_path}/aux_out were successfully created.')
        except FileExistsError: pass

        self.run_file = '#!/bin/bash\nMOPAC_DIR="/path/to/MOPAC"\nexport MOPAC_LICENSE=$MOPAC_DIR\nexport LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MOPAC_DIR\n'
        self.last_line = f'\n$MOPAC_DIR/MOPAC2016.exe ' 
        self.inner_region_size = inner_region_size
        self.number_of_snapshots = 0

    def find_atom_type(self, atom: str) -> int:
        try:
            try:
                atomic_number = atomic_properties.ATOMIC_NUMS[atom[0:2]]
                return atomic_number
            except KeyError:
                atomic_number = atomic_properties.ATOMIC_NUMS[atom[0]]
                return atomic_number
        except KeyError: raise WrongAtomTypeError(f'This atoms type is not known {atom}')

    def get_positions_in_A(self, pos: list) -> list:
        positions = []
        for p in pos:
            p = float(p) * 10
            positions.append(p)
        return positions
    
    def create_run_files(self) -> list:
        run_files_path = os.path.join(self.mopac_path, 'run_files')
        mop_files = sorted(glob(os.path.join(self.mopac_path, '*.mop')))#sorted([x for x in os.listdir(self.mopac_path) if os.path.isfile(os.path.join(self.mopac_path, x))])
        paths = []
        for mop_file in mop_files:
            name = mop_file.split('/')[-1].split('.')[0] + '.run'
            path = os.path.join(run_files_path, name)
            paths.append(path)
            with open(path,  mode='w', encoding='utf8') as w_file:
                os.system(f'chmod u+x {path}')
                w_file.write(self.run_file)
                last_line = self.last_line + os.path.join(self.mopac_path, mop_file)
                w_file.write(last_line)
        return paths
    
    def input_to_mop(self, cnf_file: str) -> None:
        '''
        The function converts gromos cnf file 
        into xyz input for MOPAC (.mop file).

        Input: cnf file
        Output: None
        '''
        name_pattern = r'_(\d+)\.'
        pos_pattern = re.compile(r'''
        ^\s+\d+\s+\w+\s+
        (\w+)\s+\d+\s+
        (\-?\d+\.\d+)\s+
        (\-?\d+\.\d+)\s+
        (\-?\d+\.\d+)\s+
        ''', re.VERBOSE)
        with open(cnf_file, 'r', encoding = 'utf8') as r_file:
            name = 'complex_snapshot_' + re.search(name_pattern, cnf_file).groups()[0] + '.mop'
            with open(os.path.join(self.mopac_path, name), 'w', encoding = 'utf8') as w_file:
                w_file.write(self.key_words)
                w_file.write(f'Generated from: {cnf_file}\n\n')
                line = r_file.readline()
                while line:
                    if line.startswith('POSITION'):
                        c = True
                        while c == True:
                            line = r_file.readline()
                            if line.startswith('END'):
                                c = False
                            else:
                                try:
                                    match = re.search(pos_pattern, line).groups()
                                    flag = 1
                                    atom = self.find_atom_type(match[0])
                                    positions = self.get_positions_in_A(match[1:])
                                    w_file.write(f'{atom:>2}{positions[0]:14.9f} {flag} {positions[1]:14.9f} {flag} {positions[2]:14.9f} {flag}\n')
                                except: pass
                    line = r_file.readline()
                #logger.info(f'File {name} was successfully written, directory: {self.mopac_path}')
        return None
    
    def _run(self, path: str) -> None:
            res = sbp.run([path], stderr=sbp.PIPE, text=True)
            if res.returncode != 0:
                raise MopacCalculationError(res.stderr)
            return None
     
    def run(self, paths: list=None)->None:
        with mp.Pool(mp.cpu_count()) as p:
            p.map(self._run, paths)
        return None
    
    def clean_up_min(self, cleaned_dir: str, output_dir: str):
        '''
        Clean up after the 1SCF the calculations after minimization.
        Mopac .arc and .out files are removed.
        '''
        #logger.info(f'Cleaning {cleaned_dir} directory...')
        os.system(f'mv {cleaned_dir}/*.aux {output_dir}')
        rmtree(cleaned_dir)

    def clean_up(self, cleaned_dir: str, output_dir: str):
        '''
        Clean up after all the calculations.
        Mopac .arc and .out files are removed.
        '''
        #logger.info(f'Cleaning {cleaned_dir} directory...')
        os.system(f'mv {cleaned_dir}/*.aux {output_dir}')
        os.system(f'rm {cleaned_dir}/*.arc')
        os.system(f'rm {cleaned_dir}/*.out')
        os.system(f'rm {cleaned_dir}/*.mop')
        rmtree(os.path.join(cleaned_dir, 'run_files'))
    
    def run_calculation(self) -> None:
        '''
        Parallelization of individual MOPAC calculations
        Clean up: .out and .arc files are removed, 
        .aux files are moved to the aux_out directory.
        '''
        if self.coord_files[0].endswith('.aux'):
            indices = []
            prev = 0
            for file in self.coord_files:
                ind = self.coord_files.index(file)
                if ind != 0:
                        self.step = self._calculate_difference(self.coord_files[ind], self.coord_files[ind - 1])
                if self.topo_path == None:
                    self.input_to_mop_buffer_freezed(file)
                else:
                    self.mopac_path, self.number_of_snapshots = ReShake(file, self.mopac_path, 
                                                                        self.number_of_snapshots, 
                                                                        self.inner_region_size, 
                                                                        self.topo_path, self.step, self.bin,
                                                                        self.imd_example, self.mk_example,
                                                                        self.lib).reshake()
                    try:
                        os.system('rm reshaken_snapshot_*')
                    except FileNotFoundError: pass
                    reshaken_cnf_files = sorted(glob(os.path.join(self.mopac_path, '*.cnf')))
                    for reshaken_cnf_file in reshaken_cnf_files:
                        self.input_to_mop(reshaken_cnf_file)
                run_paths = self.create_run_files()
                if indices == []:
                    indices.append(0)
                indices.append(len(run_paths) + prev - 1)
                prev += len(run_paths)
                self.run(run_paths)

                local_mop_path_buffer = self.input_buffer()
                run_paths_buffer = self.create_run_files_buffer_minimization(local_mop_path_buffer)
                self.run(run_paths_buffer)

                self.clean_up_min(self.mopac_path, output_dir=os.path.join(self.mopac_dir, 'minimization', 'buffer_pls_inner', 'aux_out'))
                self.clean_up_min(local_mop_path_buffer, output_dir=os.path.join(self.mopac_dir, 'minimization', 'buffer', 'aux_out'))
            np.save(os.path.join(self.mopac_dir, 'minimization', 'indices.npy'), np.array(indices, dtype=int))
            #logger.info(f'Calculations finished!')
        else:
            # inner plus buffer
            for file in self.coord_files:
                self.input_to_mop(file)
            run_paths = self.create_run_files()
            self.run(run_paths)

            # buffer
            if os.path.isdir(os.path.join(self.mopac_dir, 'buffer')):
                self.mopac_path = self.input_buffer()
                run_paths = self.create_run_files()
                self.run(run_paths)

                self.clean_up(self.mopac_path, output_dir=os.path.join(self.mopac_path, 'aux_out'))
                #logger.info(f'Calculations finished!')
            if self.mopac_path.count('buffer_pls_inner') == 1:
                self.clean_up(self.mopac_path, output_dir=os.path.join(self.mopac_path, 'aux_out'))
        return None

class Mopac_1scf_Calculation_Gromos_In(Mopac_Calculation):
    '''
    Class for running 1scf Mopac calculation (without energy minimization).
    Methods: input_to_mop
             input_buffer 
    '''
    def __init__(self, coord_files: list, inner_region_size)->None:
        if os.path.isdir('./MOPAC_results'):
            os.system('rm -r ./MOPAC_results')
        super().__init__(coord_files=coord_files, inner_region_size=inner_region_size)
        self.key_words = 'PM7 GRAD AUX(PRECISION = 9, XP, XS, XW) PRECISE 1SCF CHARGE=0\n'
        return None

    def input_buffer(self) -> str:
        '''
        The method reads the .mop files from ./MOPAC_results/buffer_pls_inner
        and creates the .mop files with the extracted buffer region in
        the ./MOPAC_results/buffer directory.

        output: path to .mop files of the buffer region: ./MOPAC_results/buffer
        '''
        inner_pls_buffer_mop_files = sorted([x for x in os.listdir(self.mopac_path) if x.endswith('.mop')])
        buffer_path = os.path.join(self.mopac_dir, 'buffer')
        try:
            os.mkdir(os.path.join(buffer_path, 'aux_out'))
            os.mkdir(os.path.join(buffer_path, 'run_files'))
        except FileExistsError: pass

        for file in inner_pls_buffer_mop_files:
            with open(os.path.join(self.mopac_path, file), mode='r', encoding='utf8') as r_file:
                num = file.split('_')[-1].split('.')[0]
                lines = r_file.readlines()
                with open(os.path.join(buffer_path, f'buffer_snapshot_{num:0>5}.mop'), mode='w', encoding='utf8') as w_file:
                    for i in range(len(lines)):
                        if i < 3 or i > self.inner_region_size + 2:
                            w_file.write(lines[i])
        self.clean_up(self.mopac_path, output_dir=os.path.join(self.mopac_path, 'aux_out'))
        return buffer_path

class Mopac_Minimization_Calculation(Mopac_Calculation):
    '''
    Class for running Mopac calculation with energy minimization
    Additional arguments: 
    freeze_surroundings: if True optimization flags for the buffer region are
                         set to 0. 
    '''
    def __init__(self, coord_file: str, freeze_buffer: bool = False):
        if os.path.isdir('./MOPAC_results'):
            os.system('rm -r ./MOPAC_results')
        super().__init__(coord_files=coord_file)
        self.key_words = 'PM7 GRAD AUX(PRECISION = 9, XP, XS, XW) CHARGE=0\n'
        self.freeze_buffer = freeze_buffer
        os.system(f'rm -r {os.path.join(self.mopac_dir, "buffer")}')

    def input_to_mop(self, cnf_file: str) -> None:
        '''
        The function converts gromos cnf file 
        into xyz input for MOPAC (.mop file).

        Input: cnf file
        Output: None
        '''
        pattern = re.compile(r'''
        ^\s+\d+\s+(\w+)\s+
        (\w+)\s+\d+\s+
        (\-?\d+\.\d+)\s+
        (\-?\d+\.\d+)\s+
        (\-?\d+\.\d+)\s+
        ''', re.VERBOSE)
        with open(cnf_file, 'r', encoding = 'utf8') as r_file:
            name = cnf_file.split('/')[-1].split('.')[0] + '.mop'
            with open(os.path.join(self.mopac_path, name), 'w', encoding = 'utf8') as w_file:
                w_file.write(self.key_words)
                w_file.write(f'Generated from: {cnf_file}\n\n')
                for line in r_file:
                    try:
                        match = re.search(pattern, line).groups()
                        mol_type = match[0]
                        flag = 1
                        if self.freeze_buffer:
                            if mol_type == 'SOLV':
                                flag = 0
                        atom = self.find_atom_type(match[1])
                        positions = self.get_positions_in_A(match[2:])
                        w_file.write(f'{atom:>2}{positions[0]:14.9f} {flag} {positions[1]:14.9f} {flag} {positions[2]:14.9f} {flag}\n')
                    except: pass
                #logger.info(f'File {name} was successfully written, directory: {self.mopac_path}')
        return None

class Mopac_1scf_Calculation_Aux_In(Mopac_Calculation):
    '''
    Class for running 1scf Mopac calculation on individual
    steps of the individual minimization calculations.
    Creates additional file system in the directory minimization
    created by class Mopac_Calculation.
    
    methods: input_to_mop
             input_buffer
             create_run_files_buffer
    '''
    def __init__(self, coord_files: list, 
                 inner_region_size: int, 
                 topo_path: str = None,
                 gromos_bin: str = None,
                 imd_example: str = None,
                 mk_example: str = None,
                 mk_lib: str = None):
        if os.path.isdir('./MOPAC_results/minimization'):
            os.system('rm -r ./MOPAC_results/minimization')
        super().__init__(coord_files=coord_files, inner_region_size=inner_region_size, 
                         topo_path=topo_path, gromos_bin=gromos_bin, imd_example=imd_example, mk_example=mk_example,
                         mk_lib=mk_lib)
        self.key_words = 'PM7 GRAD AUX(PRECISION = 9, XP, XS, XW) PRECISE 1SCF CHARGE=0\n'

    def _calculate_difference(self, acc_aux_file: str, prev_aux_file: str) -> int:
        acc = int(re.search(r'(\d+).aux', acc_aux_file).groups()[0])
        prev = int(re.search(r'(\d+).aux', prev_aux_file).groups()[0])
        return acc - prev

    def input_to_mop_buffer_freezed(self, aux_file: str) -> str:
        '''
        The method extracts coordinates of individual minimization
        steps from the given AUX file.
        '''
        file_num = int(re.search(r'(\d+).aux', aux_file).groups()[0])
        if self.mopac_path.count('snapshot') == 0:
            self.mopac_path = os.path.join(self.mopac_path, f'snapshot_{file_num:0>5}')
        else:
            pattern = f'/snapshot_{file_num-self.step:0>5}'
            self.mopac_path = self.mopac_path.replace(pattern, '')
            self.mopac_path = os.path.join(self.mopac_path, f'snapshot_{file_num:0>5}')
        try:
            os.mkdir(self.mopac_path)
        except FileExistsError: pass
       
        try:
            os.mkdir(f'{self.mopac_path}/run_files')
        except FileExistsError: pass

        with open(aux_file, mode='r', encoding='utf8') as r_file:
            line = r_file.readline()
            atoms = []
            while line:
                # ATOMS
                if line.startswith(' ATOM_EL'):
                    c = True
                    while c == True:
                        line = r_file.readline()
                        if line.startswith(' ATOM_CORE'):
                            c = False
                        else:
                            line = re.sub(' +', ' ', line).strip().split(' ')
                            for a in line:
                                atoms.append(atomic_properties.ATOMIC_NUMS[a])
                # POSITIONS 
                if line.startswith(' ATOM_X_UPDATED'):
                    with open(os.path.join(self.mopac_path, f'complex_snapshot_{self.number_of_snapshots:0>7}.mop'), mode='w', encoding='utf8') as w_file:
                        w_file.write(self.key_words)
                        w_file.write(f'Generated from: {aux_file}\n\n')
                        c = True
                        j = 0
                        while c == True:
                            line = r_file.readline()
                            if line.startswith(' HEAT_OF_FORM_UPDATED') or line.startswith(' ###########'):
                                c = False
                            else:
                                flag = 1
                                line = re.sub(' +', ' ', line.strip()).split(' ')
                                w_file.write(f'{atoms[j]}{float(line[0]):14.9f} {flag} {float(line[1]):14.9f} {flag} {float(line[2]):14.9f} {flag}\n')
                                #w_file.write(line + '\n')
                                j += 1
                    
                    self.number_of_snapshots += 1
                line = r_file.readline()
        return None #self.mopac_path
    
    def input_buffer(self) -> str:
        '''
        The method reads the .mop files from 
        ./MOPAC_results/buffer_pls_inner/minimization/buffer_pls_inner/directory_of_the_given_snapshot
        and creates the .mop files with the extracted buffer region in
        the ./MOPAC_results/buffer directory.

        output: path to .mop files of the buffer region: ./MOPAC_results/minimization/buffer/directory_of_the_given_snapshot
        '''
        inner_pls_buffer_mop_files = sorted([x for x in os.listdir(self.mopac_path) if x.endswith('.mop')])
        snapshot = self.mopac_path.split('/')[-1]
        buffer_path = os.path.join(self.mopac_dir, 'minimization', 'buffer')
        buffer_path_acc = os.path.join(buffer_path, snapshot)
        try:
            os.mkdir(os.path.join(buffer_path, 'aux_out'))
        except FileExistsError: pass
        try:
            os.mkdir(buffer_path_acc)
        except FileExistsError: pass
        try:
            os.mkdir(os.path.join(buffer_path_acc, 'run_files'))
        except FileExistsError: pass

        for file in inner_pls_buffer_mop_files:
            with open(os.path.join(self.mopac_path, file), mode='r', encoding='utf8') as r_file:
                num = file.split('_')[-1].split('.')[0]
                lines = r_file.readlines()
                with open(os.path.join(buffer_path_acc, f'buffer_snapshot_{num:0>7}.mop'), mode='w', encoding='utf8') as w_file:
                    for i in range(len(lines)):
                        if i < 3 or i > self.inner_region_size + 2:
                            w_file.write(lines[i])
        return buffer_path_acc
    
    def create_run_files_buffer_minimization(self, mopac_path) -> list:
        run_files_path = os.path.join(mopac_path, 'run_files')
        mop_files = sorted([x for x in os.listdir(mopac_path) if os.path.isfile(mopac_path + '/' + x)])
        paths = []
        for mop_file in mop_files:
            name = mop_file.split('.')[0] + '.run'
            path = os.path.join(run_files_path, name)
            paths.append(path)
            with open(path,  mode='w', encoding='utf8') as w_file:
                os.system(f'chmod u+x {path}')
                w_file.write(self.run_file)
                last_line = self.last_line + os.path.join(mopac_path, mop_file)
                w_file.write(last_line)
        return paths
    
class System_Properties():
    '''
    Class creates object which contains system propertis of the given system.
    Available properties and the corresponding methods:
    1.  atoms -> get_atoms
    2.  number of atoms -> get_num_of_atoms
    3   number of h2o molecules -> get_num_of_h2o
    4.  energy -> get_energy
    5.  normalized energy -> normalize_energy
    6.  forces -> get_forces
    7.  positions -> get_positions
    '''
    def __init__(self, compl_aux_file: str = None, buffer_aux_file: str = None, num_of_solute_atoms: int = None, reference_energies: tuple = None):
        #rint(Prepare_For_Running.minimize)
        ### PATTERNS FOR REGEX ###
        pos_pattern = re.compile(r'''
        \s+(\-?\d+\.\d+) # white space (one or more), followed by float number
        \s+(\-?\d+\.\d+)
        \s+(\-?\d+\.\d+)
        ''', re.VERBOSE)
        ##########################
        # SHARED
        self.complex_and_buffer_files: list[str] = [compl_aux_file, buffer_aux_file] 
        self.num_of_solute_atoms = num_of_solute_atoms
        self.num_of_h2o = []
        self.spin: float = None
        self.reference_energies = reference_energies

        # inner plus buffer
        self.compl_atoms = []
        self.compl_energy_kj = []
        self.compl_positions = []
        self.compl_forces = []

        # buffer
        self.buffer_atoms = []
        self.buffer_energy_kj = []
        self.buffer_positions = []
        self.buffer_forces = []

        for file in self.complex_and_buffer_files:
            test = file.split('/')[-1]
            with open(file, 'r', encoding = 'utf8') as r_file:
                line = r_file.readline()
                while line:
                    # ATOMS
                    if line.startswith(' ATOM_EL'):
                        c = True
                        ats = []
                        while c == True:
                            line = r_file.readline()
                            if line.startswith(' ATOM_CORE'):
                                c = False
                            else:
                                a = line.split()
                                ats += a
                        for a in ats:
                            if test.startswith('complex_'):
                                self.compl_atoms.append(atomic_properties.ATOMIC_NUMS[a])

                    # POSITIONS
                    if line.startswith(' ATOM_X_OPT'):
                        c = True
                        while c == True:
                            line = r_file.readline()
                            if line.startswith(' ATOM_CHARGES'):
                                c = False
                            else:
                                atom_position = np.array(re.search(pos_pattern, line).groups(), dtype=float)
                                if test.startswith('complex_'):
                                    self.compl_positions.append(atom_position)
                                else:
                                    self.buffer_positions.append(atom_position)
                    
                    # SPIN
                    if line.startswith(' TOTAL_SPIN'):
                        self.spin = float(line.split('=')[1].rstrip().replace('D', 'E'))

                    # ENERGY
                    if line.startswith(' HEAT_OF_FORMATION'):
                        energy = float(line.split('=')[1].rstrip().replace('D', 'E')) #* KCAL_TO_KJ
                        if test.startswith('complex_'):
                            self.compl_energy_kj.append(energy)
                        else:
                            self.buffer_energy_kj.append(energy)

                    # FORCES
                    if line.startswith(' GRADIENTS'):
                        c = True
                        forces = []
                        while c == True:
                            line = r_file.readline()
                            if line.startswith(' MOLECULAR_ORBITAL_OCCUPANCIES'):
                                c = False
                            else:
                                line = line.split()
                                forces += line
                        for i in range(3, len(forces) + 1, 3):
                            atom_forces = np.array(forces[i-3:i], dtype=float)
                            atom_forces = atom_forces #* KCAL_TO_KJ
                            if test.startswith('complex_'):
                                self.compl_forces.append(atom_forces)
                            else:
                                self.buffer_forces.append(atom_forces)
        
                    # READ NEXT LINE
                    line = r_file.readline()

        # inner plus buffer
        self.compl_energy_kj = np.array(self.compl_energy_kj, dtype=float)
        self.compl_atoms = np.array(self.compl_atoms, dtype=int)
        self.compl_positions = np.array(self.compl_positions)
        self.compl_forces = np.array(self.compl_forces)
        self.compl_num_of_h2o = np.array([(len(self.compl_atoms) - self.num_of_solute_atoms) // 3], dtype=int)

        # buffer
        self.buffer_positions = np.array(self.buffer_positions)

        # buffer atoms
        delay = 0 # how many atoms are on the different position in comparison with inner plus buffer region
        for i in range(len(self.compl_positions)):
            if np.array_equal(self.compl_positions[i], self.buffer_positions[i-delay]) == False: # comparison of the appropriate indices
                self.buffer_atoms.append(0) # if the positions are not equal atom is not in the buffer region -> Z = 0
                delay += 1 # adjust delay
            else:
                self.buffer_atoms.append(self.compl_atoms[i]) # atom is on the same position -> Z = Z

        self.buffer_atoms = np.array(self.buffer_atoms)
        self.buffer_num_of_h2o = np.array([len(self.buffer_atoms[self.buffer_atoms>0]) // 3], dtype = int)
        self.buffer_energy_kj = np.array(self.buffer_energy_kj, dtype=float)
        self.buffer_forces = np.array(self.buffer_forces)

    def get_spin(self) -> float:
        return self.spin

    def get_num_of_solute_atoms(self) -> int:
        return self.num_of_solute_atoms

    def get_energy(self) -> np.ndarray:
        return self.compl_energy_kj, self.buffer_energy_kj

    def get_atoms(self) -> np.ndarray:
        return self.compl_atoms, self.buffer_atoms

    def get_positions(self) -> np.ndarray:
        return self.compl_positions, self.buffer_positions
    
    def get_forces(self) -> np.ndarray:
        return self.compl_forces, self.buffer_forces
    
    def get_num_of_h2o(self) -> np.ndarray:
        return self.compl_num_of_h2o, self.buffer_num_of_h2o
    
    def get_reference_energies(self) -> float: # not finished
        mop_file = './MOPAC_buffer_pls_inner/buffer_pls_inner_00000.mop'
        solute_lines = 3 + self.num_of_solute_atoms
        mop_files = ['solute.mop', 'h2o.mop']
        with open(mop_file, 'r', encoding='utf8') as r_file:
            orig_lines = r_file.readlines()
            with open(mop_files[0], 'w', encoding='utf8') as w_file:
                w_file.writelines(orig_lines[:solute_lines])
            with open(mop_files[1], 'w', encoding='utf8') as w_file:
                w_file.writelines(orig_lines[:3])
                w_file.writelines(orig_lines[solute_lines:solute_lines+3])
    
    def get_burnn_energy(self) -> np.ndarray:
        solute_energy_vac: float = self.reference_energies[0] #* KCAL_TO_KJ # kJ/mol
        h2o_energy_vac: float = self.reference_energies[1] #* KCAL_TO_KJ # kJ/mol
        energy_comp_norm = self.compl_energy_kj - ((self.compl_num_of_h2o * h2o_energy_vac) + solute_energy_vac)
        energy_buffer_norm = self.buffer_energy_kj - (self.buffer_num_of_h2o * h2o_energy_vac)
        burnn_energy = energy_comp_norm - energy_buffer_norm
        return burnn_energy

    def get_burnn_forces(self) -> np.ndarray:
        burnn_forces = []
        delay = 0
        for i in range(len(self.compl_forces)): # iterate over all atoms within inner plus buffer region
            if self.compl_atoms[i] != self.buffer_atoms[i]: # atoms on the same indexes are not equal -> atom is not in the buffer region -> Z = 0
                burnn_forces.append(self.compl_forces[i]) # add the forces from inner plus buffer region into the appropriate index in burnn forces
                delay += 1 # adjust delay
            else: # atoms is in the buffer region
                burnn_forces.append(self.compl_forces[i] - self.buffer_forces[i - delay]) # calculate burnn forces and add them to the list
        return np.array(burnn_forces)
        
