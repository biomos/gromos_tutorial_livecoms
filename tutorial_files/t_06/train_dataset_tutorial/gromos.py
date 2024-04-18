import os, re
from pathlib import Path
import subprocess as sbp
import multiprocessing as mp
import logging as log
import atomic_properties
from glob import glob
from shutil import rmtree
import json

#logger = log.getLogger(__name__)

#formatter = log.Formatter('%(asctime)s-%(module)s-%(funcName)s-%(levelname)s: %(message)s')

#file_handler = log.FileHandler('gromos.log')
#file_handler.setLevel(log.INFO)
#file_handler.setFormatter(formatter)

#logger.addHandler(file_handler)

class CNF_File():

    def __init__(self, name: Path, title: str, timestep: str=None, position: str=None) -> None:
        self.name = name
        self.title = title
        self.ts, self.t = self._parse_timestep(timestep)
        self.position = position
        return None
    
    def _parse_timestep(self, ts_block):
        try: 
            ts, t = re.search(r'\s+(\d+)\s+(\d+\.\d+)', ts_block).groups()
            return ts, t
        except AttributeError: pass
    
    def get_timestep(self):
        return int(self.ts)
    
    def _write_ts(self, file_obj):
        file_obj.write(f"TIMESTEP\n\t{self.ts}\t{self.t}\nEND\n")
    
    def write_cnf(self):
        with self.name.open(mode='w', encoding='utf8') as w_file:
            w_file.write(self.title)
            self._write_ts(w_file)
            w_file.write(self.position)
        return None


class Extract_Buffer_Pls_Inner():
    """
    The class creates an object containing the content of the Gromos filter output files.
    Then it extracts the configurations of the QM regions of individual snapshots into
    separate files.
    
    Parameters:
    -----------

    filtered_cnfs, tuple: outputs of the filter program (can be one or more)

    Methods:
    --------

    extract: creates single cnf file for every filtered snapshot

    Returns:
    --------

    None
    """

    def __init__(self, filtered_cnfs: tuple) -> None:
        self.cwd = Path('.') 
        self.complex_files = self.cwd / 'buffer_pls_inner_region_configurations'
        self.filtered_cnfs = filtered_cnfs

        try:self.complex_files.mkdir()
        except FileExistsError: 
            rmtree(self.complex_files)
            self.complex_files.mkdir()

        self.output_file = 'buffer_pls_inner'
        self.cnf_files = []
        for filtered_cnf in self.filtered_cnfs:
            with Path(filtered_cnf).open(mode='r', encoding = 'utf8') as r_file:
                self.cnf_files.append(r_file.readlines())
                #logger.info(f'The input file {input} was read in.')
        return None

    def _lines_generator(self, lines: list) -> str:
        for line in lines:
            yield line

    def extract(self) -> None:
        i = 0
        for lines in self.cnf_files:
            generator = self._lines_generator(lines) # generator of the individual lines  
            header = '' # TITLE block
            for header_line in generator:
                header += header_line
                if header_line == 'END\n': # END of the block
                    break

            line = next(generator)
            while line:
                if line.startswith('TIMESTEP'):
                    cnf_lines = [] # lines of one .cnf file
                    c = True
                    while c: # three blocks per file
                        if line.startswith('GENBOX'):
                            c = False
                        else:
                            cnf_lines.append(line)
                        line = next(generator)
                    output_file_path = self.complex_files / f'{self.output_file}_{i:0>5}.cnf'
                    with output_file_path.open(mode='w', encoding = 'utf8') as w_file:
                        w_file.write(header) # TITLE block
                        for line in cnf_lines: # writing lines of the single file
                            w_file.write(line)
                        i += 1
                        #logger.info(f'cnf file {self.output_file}_{i:0>5}.cnf was written.')
                try:
                    line = next(generator)
                except StopIteration: break
        return None
    
class Extract_Buffer_Pls_Inner_Adaptive_Sampling(Extract_Buffer_Pls_Inner):
    def __init__(self, filtered_cnfs: tuple, simulation_dirs: tuple) -> None:
        super().__init__(filtered_cnfs)
        self.simuldirs = simulation_dirs

    def _key(self, omd_file: Path)->int:
        return int(re.search(r'(\d+)\.omd', str(omd_file)).groups()[0])
    
    def find_omd_files(self, simuldir: Path):
        return sorted(simuldir.glob('*.omd'), key=self._key)
    
    def find_nnvalid_lines(self, lines: list):
        pattern = re.compile(r'step\s+(\d+)\s+:')
        steps = []
        for line in lines:
            if line.startswith('   NOTICE NN Worker : Deviation from validation model above threshold in step'):
                steps.append(int(re.search(pattern, line).groups()[0]) + 1)
        return steps
                
    def find_steps_above_nnvalid_thr(self)->None:
        self.steps_above_nnvalid_thr = []
        for simuldir in self.simuldirs:
            simuldir = Path(simuldir)
            omd_files = self.find_omd_files(simuldir=simuldir)
            for omd_file in omd_files:
                with omd_file.open(mode='r', encoding='utf8') as r_file:
                    self.steps_above_nnvalid_thr += self.find_nnvalid_lines(r_file.readlines())
        return None
    
    def _extract_title(self, filter_cnf_generator)->str:
        title = '' # TITLE block
        for title_line in filter_cnf_generator:
            title += title_line
            if title_line == 'END\n': # END of the block
                break
        return title
    
    def _extract_property(self, line, filter_cnf_generator, property:str):
        property_str = '' 
        while line:
            if line.startswith(property.upper()):
                while True:
                    property_str += line
                    if line.startswith('END'):
                        break
                    line = next(filter_cnf_generator)
            if property_str != '':
                return property_str
            else: 
                try: line = next(filter_cnf_generator)
                except StopIteration: break

    
    def extract_individual_cnfs(self):
        c = 0
        for filtered_cnf in self.filtered_cnfs:
            with Path(filtered_cnf).open(mode='r', encoding='utf8') as r_file:
                title = self._extract_title(r_file)
                line = next(r_file)
                while line:
                    time_step = self._extract_property(line, r_file, property='TIMESTEP')
                    try: line = next(r_file)
                    except StopIteration: break
                    configuration = self._extract_property(line, r_file, property='POSITION')
                    output_file_path = self.complex_files / f'{self.output_file}_{c:0>5}.cnf'
                    cnf_file = CNF_File(name=output_file_path, title=title, timestep=time_step, position=configuration)
                    ts = cnf_file.get_timestep()
                    if ts == self.steps_above_nnvalid_thr[c]:
                        cnf_file.write_cnf()
                        c += 1
                    line = next(r_file)
    
class Extract_Buffer_Pls_Inner_Time_Clustering(Extract_Buffer_Pls_Inner):
    """
    The class creates an object containing the content of the Gromos filter output files.
    Then it extracts the configurations of the QM regions of individual snapshots into
    separate files.

    Parameters:
    -----------

    filtered_cnfs, tuple: outputs of the filter program (can be one or more)
    simulation_dirs, list: simulation directory (directory where to search for omd files)
    save_json, bool: if True, control json file with time clusters is saved

    Methods:
    --------

    time_clustering: selects snapshots based on the time occurence (if there are more snapshots in 
    the row only one with the highest NN valid deviation is selected)

    extract: creates single cnf file for every filtered snapshot if it is selected by time clustering

    Returns:
    --------

    None
    """

    def __init__(self, filtered_cnfs: tuple, simulation_dirs: list, save_json: bool=False)->None:
        super().__init__(filtered_cnfs)
        self.sim_dirs = simulation_dirs
        self.save_json = save_json
        
    def _time_clustering(self, omd_files: list):
        pattern = re.compile(r'''
        \s+WARNING\s+NN\s+Worker\s+\:\s+Deviation\s+from\s+validation\s+model\s+above\s+threshold\s+in\s+step\s+
        (\d+)\s+\:\s+(\-?\d+\.\d+)
        ''', re.VERBOSE)
        clusters = {}
        time_steps = []
        for i in range(len(omd_files)):
            key = re.search(r'(\d+).omd', omd_files[i]).groups()[0]
            with open(omd_files[i], 'r', encoding='utf8') as r_file:
                t_steps = [] # time steps
                devs = [] # NN valid deviations
                prev_dev = 0
                prev_step = 0
                for line in r_file:
                    try:
                        match = re.search(pattern, line).groups()
                        dev = float(match[1]) # NN valid deviation
                        step = int(match[0]) + 1 # time step
                        if prev_dev < abs(dev) and step == prev_step + 1:
                            if len(devs) != 0:
                                t_steps.pop(-1)
                                time_steps.pop(-1)
                                devs.pop(-1)
                            t_steps.append(step)
                            time_steps.append(step)
                            devs.append(dev)
                            prev_dev = abs(dev)
                            prev_step = step
                        elif step > prev_step + 1:
                            t_steps.append(step)
                            time_steps.append(step)
                            devs.append(dev)
                            prev_dev = abs(dev)
                            prev_step = step
                        else:
                            prev_step = step
                    except AttributeError: pass
                if len(t_steps) != 0:
                    clusters[key] = (tuple(t_steps), tuple(devs))
                 
        return time_steps, clusters
    
    def time_clustering(self)->None:
        self.time_clusters = {}
        control_dict = {}
        for i in range(len(self.sim_dirs)):
            key = i + 1
            sim_path = os.path.relpath(self.sim_dirs[i])
            omd_files = glob(os.path.join(sim_path, '*.omd'))
            t_steps, clusters = self._time_clustering(omd_files)
            self.time_clusters[key] = t_steps
            control_dict[key] = clusters
        if self.save_json:
            with open('time_clusters_from_as_burnn_simulation.json', 'w', encoding='utf8') as w_file:
                json.dump(control_dict, w_file)
        return None
    
    def _determine_timestep(self, line: str) -> int:
        ts = int(re.sub(' +', ' ', line).strip().split(' ')[0])
        return ts

    def extract(self) -> None:
        i = 0
        for idx in range(len(self.cnf_files)):
            required_timesteps = self.time_clusters[idx+1]
            generator = self._lines_generator(self.cnf_files[idx]) # generator of the individual lines  
            header = '' # TITLE block
            for header_line in generator:
                header += header_line
                if header_line == 'END\n': # END of the block
                    break

            line = next(generator)
            j = 0
            while line:
                if line.startswith('TIMESTEP'):
                    cnf_lines = [] # lines of one .cnf file
                    c = True
                    while c: # three blocks per file
                        if line.startswith('GENBOX'):
                            c = False
                        else:
                            if len(cnf_lines) == 1:
                                time_step = self._determine_timestep(line)
                                try:
                                    if time_step != required_timesteps[j]:
                                        break
                                    else:
                                        j += 1
                                except IndexError: break
                            cnf_lines.append(line)
                        line = next(generator) 
                    if len(cnf_lines) > 1:
                        with open(f'./buffer_pls_inner_region_configurations/{self.output_file}_{i:0>5}.cnf', 'w', encoding = 'utf8') as w_file:
                            w_file.write(header) # TITLE block
                            for line in cnf_lines: # writing lines of the single file
                                w_file.write(line)
                            i += 1
                        #logger.info(f'cnf file {self.output_file}_{i:0>5}.cnf was written.')
                try:
                    line = next(generator)
                except StopIteration: break
        return None
    
class ReShakeError(Exception):
    pass

class ReShake():

    def __init__(self, aux_file: str=None, 
                 mopac_path: str=None, 
                 number_of_snapshots: int=None, 
                 inner_region_size: int=None, 
                 topo_path: str=None, 
                 step: int=None,
                 gromos_bin: str=None,
                 imd_example: str=None,
                 mk_example: str=None,
                 mk_lib: str=None) -> None:
        self.aux_file = aux_file
        self.mopac_path = mopac_path
        self.number_of_snapshots = number_of_snapshots
        self.inner_region_size = inner_region_size
        self.num_of_atoms = None
        self.topo_path = topo_path
        self.step = step
        self.bin = gromos_bin
        self.imd_example = imd_example
        self.mk_example = mk_example
        self.lib = mk_lib

    def aux_to_cnf(self) -> str:
        '''
        The method extracts coordinates of individual minimization
        steps from the given AUX file.
        '''
        file_num = int(re.search(r'(\d+).aux', self.aux_file).groups()[0])
        if self.mopac_path.count('snapshot') == 0:
            self.mopac_path = os.path.join(self.mopac_path, f'snapshot_{file_num:0>5}')
        else:
            pattern = f'/snapshot_{file_num-self.step:0>5}'
            self.mopac_path = self.mopac_path.replace(pattern, '')
            self.mopac_path = os.path.join(self.mopac_path, f'snapshot_{file_num:0>5}')
        
        if os.path.isdir(self.mopac_path):
            rmtree(self.mopac_path)
        try:
            os.mkdir(self.mopac_path)
        except FileExistsError: pass
        try:
            os.mkdir(f'{self.mopac_path}/run_files')
        except FileExistsError: pass

        with open(self.aux_file, mode='r', encoding='utf8') as r_file:
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
                            line = line = re.sub(' +', ' ', line).strip().split(' ')
                            for a in line:
                                atoms.append(atomic_properties.ATOMIC_NUMS[a])
                    if self.num_of_atoms == None:
                        self.num_of_atoms = len(atoms)
                # POSITIONS 
                if line.startswith(' ATOM_X_UPDATED'):
                    with open(os.path.join(self.mopac_path, f'complex_snapshot_{self.number_of_snapshots:0>7}.cnf'), mode='w', encoding='utf8') as cnf_file:
                        cnf_file.write('TITLE\n')
                        cnf_file.write(f'\tGenerated from: {self.aux_file}\n')
                        cnf_file.write('END\n')
                        cnf_file.write('POSITION\n')
                        c = True
                        j = 0
                        while c == True:
                            line = r_file.readline()
                            if line.startswith(' HEAT_OF_FORM_UPDATED') or line.startswith(' ###########'):
                                c = False
                            else:
                                line = re.sub(' +', ' ', line.strip()).split(' ')
                                cnf_file.write(f'{j+1:5}   UNL    *{j+1:8}{float(line[0])/10:16.9f}{float(line[1])/10:16.9f}{float(line[2])/10:16.9f}\n')
                                #w_file.write(line + '\n')
                                j += 1
                        cnf_file.write('END')
                    self.number_of_snapshots += 1
                line = r_file.readline()
        return self.mopac_path
    
    def _parse_example(self, file):
        with open(file, mode= 'r', encoding='utf8') as file:
            return file.readlines()
        
    def create_imd_file(self)->None:
        example_imd = self._parse_example(self.imd_example)
        with open(os.path.join(self.mopac_path, 'reshake.imd'), mode='w', encoding='utf8') as imd_file:
            for i in range(len(example_imd)):
                if re.search(r'#\s+NPM', example_imd[i]) != None:
                    example_imd[i+1] = f'         1       {(self.num_of_atoms - self.inner_region_size)//3}\n'
                elif re.search(r'#\s+NEGR\s+NRE', example_imd[i]) != None:
                    example_imd[i+1] = f'     2       {self.inner_region_size}        {self.num_of_atoms}\n'
                imd_file.write(example_imd[i])
        return None
    
    def create_arg_files(self)-> None:
        example_arg = self._parse_example(self.mk_example)
        cnf_files_names = glob(os.path.join(self.mopac_path, '*.cnf'))
        for cnf_file_name in cnf_files_names:
            cnf_file = cnf_file_name.split('/')[-1]
            arg_name = 'mk_' + cnf_file.split('.')[0].split('_')[2] + '.arg'
            with open(os.path.join(self.mopac_path, arg_name), 'w') as w_file:
                for line in example_arg:
                    if line.startswith('@sys'):
                        line = f'@sys\treshaken_snapshot_{cnf_file.split(".")[0].split("_")[2]}\n'
                        w_file.write(line)
                    elif line.startswith('@dir'):
                        line = f'@dir\t{self.mopac_path}\n'
                        w_file.write(line)
                    elif line.startswith('  topo'):
                        line = f'  topo\t{self.topo_path}\n'
                        w_file.write(line)
                    elif line.startswith('  coord'):
                        line = f'  coord\t{cnf_file}\n'
                        w_file.write(line)
                    elif line.startswith('@template'):
                        line = f'@template\t{self.lib}\n'
                        w_file.write(line)
                    elif line.startswith('@bin'):
                        line = f'@bin\t{self.bin}\n'
                        w_file.write(line)
                    else:
                        w_file.write(line)
        return None

    def create_run_files(self, mk_file):
        res = sbp.run(['mk_script', '@f', f'{mk_file}'], stderr=sbp.PIPE, text=True)
        if res.returncode != 0:
            raise ReShakeError(res.stderr)

    def run_reshake(self, run_file):
        res = sbp.run([run_file], stderr=sbp.PIPE, text=True)
        if res.returncode != 0:
            raise ReShakeError(res.stderr)

    def reshake(self):
        self.aux_to_cnf()
        self.create_imd_file()
        self.create_arg_files()
        mk_files = glob(os.path.join(self.mopac_path, '*.arg'))
        with mp.Pool(mp.cpu_count()) as p:
            p.map(self.create_run_files, mk_files)
        
        run_files = glob(os.path.join(self.mopac_path, '*.run'))
        with mp.Pool(mp.cpu_count()) as p:
            p.map(self.run_reshake, run_files)

        os.system(f'rm {self.mopac_path}/*.run')
        os.system(f'rm {self.mopac_path}/*.arg')
        os.system(f'rm {self.mopac_path}/*.imd')
        os.system(f'rm {self.mopac_path}/*.omd')
        os.system(f'rm {self.mopac_path}/mess')
        os.system(f'rm {self.mopac_path}/complex_*')

        final_cnf_files = glob(os.path.join(self.mopac_path, '*.cnf'))
        for final_cnf in final_cnf_files:
            split = final_cnf.split('/')
            name = split[-1]
            name = re.search(r'(\D+)_(\D+)_(\d+)_\d(.cnf)', name).groups()
            name = '_'.join(name[:3]) + name[-1]
            path_list = split[1:-1]
            path_list.append(name)
            path = os.path.join('/', *path_list)
            os.system(f'mv {final_cnf} {path}')

        return self.mopac_path, self.number_of_snapshots

class Create_Qmm_File():
    def __init__(self, input_file,
                 output_file_name,
                 qm_program: str = None,
                 nnmodel_path: str = None, 
                 nnmodel_pars: tuple = (0, 1),
                 nnvalid_path: str = None,
                 nnvalid_pars: tuple = (1, 5.0, 0.0),
                 nncharge_path: str = None,
                 qmzone_size: int = None,
                 qmzone_pars: tuple = (0, 1),
                 bufferzone_pars: tuple = (0, 1, 0.5)) -> None:
        
        self.unit_conversion = {'mopac': [0.1, 4.184, -41.84, 1.0]}
        self.input_file = input_file
        self.name = output_file_name
        self.qm_program = qm_program
        self.nnmodel_path = nnmodel_path
        self.nnvalid_path = nnvalid_path
        self.nncharge_path = nncharge_path
        self.qmzone_size = qmzone_size
        self.nnmodel_pars = nnmodel_pars
        self.nnvalid_pars = nnvalid_pars
        self.qmzone_pars = qmzone_pars
        self.bufferzone_pars = bufferzone_pars
        return None
    
    def read_input_file(self) -> list:
        pattern = re.compile(r'''
        \s+(\d+)\s+
        (\w+)\s+
        (\w+)\s+(\d+)\s+
        ''', re.VERBOSE)
        try:
            with open(self.input_file, mode='r', encoding='utf8') as r_file:
                qmzone = []
                bufferzone = []
                line = r_file.readline()
                while line:
                    if line.startswith('POSITION'):
                        c = True
                        i = 0
                        while c:
                            if line.startswith('END'):
                                c = False
                            try:
                                match = list(re.search(pattern, line).groups())
                                match.append(atomic_properties.ATOMIC_NUMS[match[2][0]])
                                match.append(0)
                                if i < self.qmzone_size:
                                    qmzone.append(match)
                                    i += 1
                                else:
                                    bufferzone.append(match)
                            except AttributeError: pass
                            line = r_file.readline()
                    line = r_file.readline()
                return qmzone, bufferzone

        except FileNotFoundError: print(f' A file {self.input_file} does not exists.')
    
    def create_qmm(self):
        qm_zone, buffer_zone = self.read_input_file()
        unit_conversion = self.unit_conversion[self.qm_program]
        with open(self.name + '.qmm', mode='w', encoding='utf8') as w_file:
            w_file.write('TITLE\nqmmm specification file\nEND\nQMUNIT\n')
            w_file.write('# QMULEN: conversion factor for distance units of the model\n')
            w_file.write('# QMUENE: conversion factor for energy units of the model\n')
            w_file.write('# QMUFOR: conversion factor for force units of the model\n')
            w_file.write('# QMULEN    QMUENE    QMUFOR    QMUCHR\n')
            w_file.write(f'{unit_conversion[0]:>8}{unit_conversion[1]:>8}{unit_conversion[2]:>8}{unit_conversion[3]:>8}\nEND\n')
            w_file.write('NNDEVICE\ncuda\nEND\n')
            w_file.write(f'NNMODEL\n{self.nnmodel_path}\n{self.nnmodel_pars[0]} {self.nnmodel_pars[1]}\nEND\n')
            w_file.write(f'NNVALID\n{self.nnvalid_path}\n{self.nnvalid_pars[0]} {self.nnvalid_pars[1]} {self.nnvalid_pars[2]}\nEND\n')
            if self.nncharge_path != None:
                w_file.write(f'NNCHARGE\n{self.nncharge_path}\n0\nEND\n')
            w_file.write(f'QMZONE\n{self.qmzone_pars[0]} {self.qmzone_pars[1]}\n')
            w_file.write('# RESIDUE   ATOM     QMI   QMZ   QMLI\n')
            for line in qm_zone:
                w_file.write(f'{line[0]:>4} {line[1]:<6} {line[2]:<5} {line[3]:>6}{line[4]:>6}{line[5]:>6}\n')
            w_file.write('END\n')
            w_file.write(f'BUFFERZONE\n{self.bufferzone_pars[0]} {self.bufferzone_pars[1]} {self.bufferzone_pars[2]}\n')
            for line in buffer_zone:
                w_file.write(f'{line[0]:>4} {line[1]:<6} {line[2]:<5} {line[3]:>6}{line[4]:>6}{line[5]:>6}\n')
            w_file.write('END\n')