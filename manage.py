""" Python interface to read and set-up the simulation(s) according
to the specifications given in the lib/config.yaml file.
"""
try:
    import os
except ImportError:
    print(f"[ FAIL ]...unable to import the python package 'os'.")
    sys.exit()
try:
    import sys
except ImportError:
    print(f"[ FAIL ]...unable to import the python package 'sys'.")
    sys.exit()
try:
    import yaml
except ImportError:
    print(f"[ FAIL ]...unable to import the python package 'yaml'.")
    sys.exit()
try:
    import shutil
except ImportError:
    print(f"[ FAIL ]...unable to import the python package 'shutil'.")
    sys.exit()

# --------------------------------------------------------------------------------------------------

GCC = "g++"
PCC = "pgc++"
OPT = "-fast -O3"
MCOREOPT = "-O3 -fopenmp"
STD = "-std=c++0x"

# --------------------------------------------------------------------------------------------------

# Job Compilation options specifier decalrations
ACC_OPT = '--acc'
PROFILE_OPT = '--prof'
MULTI_CORE_OPT = '--mcore'
SINGLE_CORE_OPT = '--score'

# --------------------------------------------------------------------------------------------------

#Job submission options
LOC = "loc"
CONDOR = "condor"
SLURM = "slurm" 

# --------------------------------------------------------------------------------------------------

presets_filename = "presets.hpp"

# --------------------------------------------------------------------------------------------------

TARGET = "main"
OBJECT = "main.o"
SOURCE = "main.cpp"

# --------------------------------------------------------------------------------------------------


file_exists_warning_prompt = ("[ WARNING ]...The existing directrories and will be removed. " + 
                            "Do you want to continue? (y/n): "
                              )
# --------------------------------------------------------------------------------------------------

instructions = ("\n" + 
    f"""                              HELP
    Run :
    -----
    +---------------------------------------------------------------+
    | $python manage.py  [compile_opt] [scheme] (--s [submit_opt])  |
    +---------------------------------------------------------------+
    Note: ( ... ) is optional

    [ compile_opt ] :
    -----------------
        {SINGLE_CORE_OPT.ljust(10)}: Single core job (assumes {GCC} is installed).
        {MULTI_CORE_OPT.ljust(10)}: Multicore job (assume {GCC} and OpenMP are installed).
        {ACC_OPT.ljust(10)}: GPU accelerated job (by default assumes {PCC} compiler is installed. 
        {' '*11} Else, edit the PCC variable in the manage.py to point to the correct compiler with 
        {' '*11} OpenACC support).

    [submit_opt] :
    --------------
        {LOC}    : Local submission
        {CONDOR} : Submit with condor
        {SLURM}  : Submit with slurm 
    eg:
        Run FV scheme with gpu & submit with slurm : $ python manage.py  {ACC_OPT} fv --s {SLURM}
        Run FD scheme with gpu & submit locally    : $ python manage.py  {ACC_OPT} fd --s {LOC}
    """
    )
# --------------------------------------------------------------------------------------------------

proj_dir = os.getcwd()
resources_dir = os.path.join(proj_dir, "src")
presets_file = os.path.join(resources_dir, presets_filename)

# Configuration file for the collection of jobs.
batch_configs_file = os.path.join(proj_dir, "cosenu_configs.yaml")

#List of the jobs will be stored in this txt file.
jobs_list_file = "jobs_list.txt"

# A job.config file will be added to each job folder
job_config_file = ""

config_file = "job.config"

# For submitting with condor
condor_submission_file = "submit.jdl"


# --------------------------------------------------------------------------------------------------

def write_dict(dct, path):
    keys = dct.keys()
    maxlen = max([len(key) for key in keys])
    with open(path, 'w') as f:
        for key in keys:
            f.write(f'{str(key).ljust(maxlen+1)} : {str(dct[key])} \n')

# --------------------------------------------------------------------------------------------------


def export_job_configs(config_dict, path):
    with open(path, 'w') as f:
        for entry in config_dict:
            f.write(f'{entry["id"]}, {entry["file"]} \n')

# --------------------------------------------------------------------------------------------------

def rm(path):
    if os.path.exists(path):
        stat = os.remove(path)
        # print(f'{path} removed.')
    else:
        print(f"{path} does not exist.")

# --------------------------------------------------------------------------------------------------

def configure(scheme="fv"):
    """Subroutine to read the configure file and make necessary folders"""
    proj_dir = os.getcwd()
    scheme_dir_path = ""
    with open(batch_configs_file) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    nzs = config["nzs"]
    nvzs = config["nvzs"]
    CFLS = config["CFLS"]
    z = config["zrange"]
    v0 = config["v0"]
    v1 = config["v1"]
    end_time = config["end_time"]
    nanalyze = config["n_analyze"]
    scheme_dir = os.path.join(proj_dir,  config[f'folder_{scheme}'])
    boundary = config['boundary']

    with open(presets_file, "w") as presets:
        presets.write("#if not defined(__PRESETS__)\n")
        presets.write("#define __PRESETS__\n\n")
        if boundary == 'open':
            presets.write("#define OPEN_BC\n")
        else:
            presets.write("#define PERIODIC_BC\n")

        if (scheme == "fd"):
            presets.write("#define FD\n")
            presets.write("#define KO_ORD_3\n")
            presets.write("#define ADVEC_CENTER_FD\n")
            
        if (scheme == "fv"):
            presets.write("#define FV\n")

        if config["vac_osc_on"] == True:
            presets.write("#define VAC_OSC_ON\n")
        if config["mat_osc_on"]==True:
            presets.write("#define MAT_OSC_ON\n")
        if config["collective_osc_on"] == True:
            presets.write("#define COLL_OSC_ON\n")
        if config["advection_off"] == True:
            presets.write("#define ADVEC_OFF\n")
        presets.write("\n")

        presets.write("#endif // __PRESETS__\n")
        
    # Creating necessary folders and copying files
    config_list = []

    continue_if_exist = "N"

    try:
        os.mkdir(scheme_dir)
        print(f"[ OK ]...Creating {scheme_dir}")
    except FileExistsError:
        print(f"[ EXISTS ]...{scheme_dir}")
        flag = True
        while flag:
            continue_if_exist = input(file_exists_warning_prompt).upper()
            if continue_if_exist == "N":
                print(f"[ OK ]...exiting the {scheme} setup.")
                sys.exit()
            elif continue_if_exist == "Y":
                print(f"[ OK ]...Setting up the projects for {scheme}.")
                shutil.rmtree(scheme_dir)
                os.mkdir(scheme_dir)
                flag = False
            else:
                print(f"[ FAIL ]...Invalid input {continue_if_exist}\n Re-enter.")
                

    if (scheme == 'fd'):
        gz = 2
    elif (scheme == 'fv'):
        gz = 4
        
    # Loding config values from config file
    for i, nz in enumerate(nzs):
        for nvz in nvzs:
            for CFL in CFLS:

                z0 = z[0]
                z1 = z[1]
                dz = (z1-z0)/nz
                dt = abs(CFL*dz/v1)
                N_ITER = int(end_time/dt)

                if nanalyze != 0:
                    ANAL_EVERY = int(N_ITER/nanalyze)
                    if (ANAL_EVERY < 1):
                        ANAL_EVERY = 1
                else:
                    ANAL_EVERY = N_ITER + 1

                current_config = {
                    'scheme': scheme.upper(),
                    'z0': z[0],
                    'z1': z[1],
                    'nvz': nvz,
                    'CFL': CFL,
                    'dz' : dz,
                    'dt' : dt,
                    'gz': gz,
                    'v0': v0,
                    'v1': v1,
                    'nz': nz,
                    'pmo': config['pmo'],
                    'omega': config['omega'],
                    'theta': config['theta'],
                    'mu': config['mu'],
                    'perturbation_size' : config['perturbation_size'],
                    'N_ITER': N_ITER,
                    'ANAL_EVERY': ANAL_EVERY,
                    'n_dump_rho' : config['n_dump_rho'],
                }

                # ID for the job
                config_id = f"{nz}_{nvz}_{CFL}"

                # Making directory with name ID
                model_dir = os.path.join(scheme_dir, config_id)

                # Config file for the job "ID"

                config_list.append({"id": config_id, "file": config_file})

                try:
                    dir = os.mkdir(model_dir)
                    print(f"[ OK ]...Created {model_dir}")
                except FileExistsError:
                    print(f"[ EXISTS ]...{model_dir}")

                config_path = os.path.join(model_dir, config_file)
                write_dict(current_config, config_path)

    job_config_file = os.path.join(scheme_dir, jobs_list_file)
    export_job_configs(config_list, job_config_file)

    os.chdir(proj_dir)

    return "success", scheme_dir
# ---------------------------------------------------------------------------------------------------
def compi(comp_opt = "--acc"):
    """To compile the code"""

    with open(batch_configs_file, "r") as con:
        configs = yaml.load(con, Loader=yaml.FullLoader)
    modules = configs["modules"]

    pwd = os.getcwd()
    os.chdir(resources_dir)
    exe_path = os.path.join(os.getcwd(), TARGET)

    rm(os.path.join(os.getcwd(), OBJECT))
    rm(os.path.join(os.getcwd(), TARGET))
    comp_stat = None

    if comp_opt == "--score":
        comp_stat = os.system(f"{GCC} {STD} -o {TARGET} {SOURCE}")
        if comp_stat == 0:
            print(f"{TARGET} generated for singlecore job")
        else:
            return None, None
    elif comp_opt == "--mcore":
        comp_stat = os.system(f"{GCC} {STD} {MCOREOPT} -o {TARGET} {SOURCE}")
        if comp_stat == 0:
            print(f"{TARGET} generated for muti-core job")
        else:
            return None, None
    elif comp_opt == "--acc":
        #print('dsdfsfsdf')
        with open ("compile.sh", "w") as f:
            f.write("#!/bin/bash" + "\n")
            f.write("module purge" + "\n")
            for module in modules:
                f.write(f"module load {module} \n")
            f.write("make" + "\n")
        os.system("chmod +x compile.sh")
        with open ("Makefile", "w") as f:
            f.write("#!/bin/bash" + "\n")
            f.write("all:" + "\n")
            f.write("\t" + f"{PCC} {OPT} -acc -Minfo=accel -ta=tesla:managed -o {TARGET} {SOURCE}" + "\n")

        comp_stat = os.system(f"./compile.sh")

        if comp_stat == 0:
            print(f"{TARGET} generated for accelerated job")
        else:
            return None, None
    else:
        print(f"Invalid compilation option {comp_opt}")
    
    os.chdir(pwd)
    stat = "failed"
    if comp_stat == 0:
        stat = "success"
    return stat, exe_path
# --------------------------------------------------------------------------------------------------

def cp_exe(scheme_dir_path, exe_path):
    """Copy executable from the compilation to the folders created(by configure)."""
    jobs_list = []
    with open(os.path.join(scheme_dir_path, jobs_list_file), 'r') as f:
        for line in f.readlines():
            job = line.strip("\n").split(",")[0]
            jobs_list.append(line.strip("\n").split(",")[0])
            dst = os.path.join(scheme_dir_path, job, TARGET)
            if os.path.exists(dst):
                os.system(f"rm {dst}")
                print(f"removed {dst}")
                shutil.copy(exe_path, dst)
            else:
                shutil.copy(exe_path, dst)
            print (f"copied {exe_path} to {dst}")
                
    return "success", jobs_list
    
# --------------------------------------------------------------------------------------------------

def run(jobs_list, scheme_dir_path, submit_mod):
    pwd = os.getcwd()
    for job in jobs_list:
        os.chdir(os.path.join(scheme_dir_path, job))
        if submit_mod == "loc":
            os.system(f"./{TARGET} --id {job} --conf {config_file}")

        elif submit_mod == "condor":

            print (f"Submitting {job} with {submit_mod}")

            with open (os.path.join(pwd, batch_configs_file), 'r') as f:

                configs = yaml.load(f, Loader=yaml.FullLoader)

            jdl_configs = configs["condor_requests"]

            with open (condor_submission_file, 'w') as f:
                
                f.write ("universe = vanilla" + "\n")
                f.write (f"ID = {job}" + "\n")
                f.write (f"request_cpus = {jdl_configs['ncpu']}" + "\n")
                f.write (f"request_memory = {jdl_configs['ram']}" + "\n")
                f.write (f"request_disk = {jdl_configs['storage']}" + "\n")
                f.write ("error = $(ID).err" + "\n")
                f.write ("output = $(ID).out" + "\n")
                f.write ("log = $(ID).log" + "\n")
                f.write (f"executable = {TARGET}" + "\n")
                f.write (f"arguments = --id {job} --conf {config_file}" + "\n")
                f.write ("should_transfer_files = yes" + "\n")
                f.write (f"transfer_input_files = {config_file}" + "\n")
                f.write ("queue 1")

            os.system(f"condor_submit {condor_submission_file}")

        elif submit_mod == "slurm":

            print (f"Submitting {job} with {submit_mod}")

            with open (os.path.join(pwd, batch_configs_file), 'r') as f:
                configs = yaml.load(f, Loader=yaml.FullLoader)
            
            slurm_requests = configs["slurm_requests"]
            



            with open("sjob.submit", "w") as f:
                f.write("#!/bin/bash" + "\n")
                
                f.write(f"#SBATCH --job-name={job}" + "\n")
                if slurm_requests["GPU"]["enable"]:
                    f.write(f"#SBATCH --partition={slurm_requests['GPU']['partition']}" + "\n")
                    f.write(f"#SBATCH --gres=gpu:{slurm_requests['GPU']['ngres']}" + "\n")
                
                if slurm_requests["CPU"]["enable"]:
                    f.write(f"#SBATCH --cpus-per-task={slurm_requests['CPU']['cpus-per-task']}" + "\n")
                    #f.write(f"module {slurm_requests['CPU']['purge']}" + "\n")
                    f.write(f"#SBATCH --partition={slurm_requests['CPU']['partition']}" + "\n")
                for module in configs['modules']:
                    f.write(f"module load {module}" + "\n")

                f.write(f"srun ./{TARGET} --id {job} --conf {config_file}" + "\n")
        
            os.system("sbatch sjob.submit")
        os.chdir(pwd)
        
    return "success"

# --------------------------------------------------------------------------------------------------

def main(mode, scheme, submit_mod = "loc", do_submit = False):
    conf_stat = "failed"  # Configuration status
    compi_stat = "failed" # Compilation status
    cp_stat = "failed"    # Copy status
    run_stat = "failed"   # Run status


    """Configure"""
    conf_stat, scheme_dir_path = configure(scheme)
    if conf_stat == "success":
        """Compile"""
        compi_stat, exe_path = compi(mode)
    else:
        print("Configuration failed")
        return

    if compi_stat == "success":
        """Copy executable to the folders"""
        cp_stat, jobs_list = cp_exe(scheme_dir_path, exe_path)
    else:
        print("Compilation failed")
        return
    
#Uncomment this section to run the jobs automatically.
    if (cp_stat == "success"):
        if do_submit:
            """Run jobs"""
            run_stat = run(jobs_list, scheme_dir_path, submit_mod)
        else:
            print ("Did not ask to submit the job.")
            return 
    else:
        if (not do_submit):
            print(f"Copying executable failed.")
        else:
            print ("You chose not to submit the job.")
        return
        
    if run_stat == "success":
        print("SUCCESS")
    else:
        print("FAILED")

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    PWD = os.getcwd()
    scheme = None
    mode = None
    is_scheme = False
    is_mode = False
    submit_mod = "loc"
    do_submit = False
    
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv)):

            if sys.argv[i] == "--help":
                print(instructions)
                exit(0)

            elif (sys.argv[i] == SINGLE_CORE_OPT) or (sys.argv[i] == MULTI_CORE_OPT) or (sys.argv[i] == ACC_OPT):
                mode  = sys.argv[i]                
                i += 1
                is_mode = True
                print (f"compilation mode set to {mode}")

            elif sys.argv[i] == 'fv' or sys.argv[i] == 'fd':
                scheme = sys.argv[i]
                i += 1
                is_scheme = True
                print(f"Simulation scheme is set to {scheme}")

            elif sys.argv[i] == "--s":
                submit_mod = sys.argv[i+1]
                do_submit = True
                print(f"Submission option set to {submit_mod}")
            

    if not (is_mode and is_scheme):
        print ("You have not specified either the simulation [scheme] or compilation [opt]")
        print("Run manage.py --help for instructions.")
        sys.exit()

    if (not is_scheme and not is_mode):
        print(instructions)
    else:
        main(mode, scheme, submit_mod = submit_mod, do_submit= do_submit)
        os.chdir(PWD)
