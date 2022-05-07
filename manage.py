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

instructions = """ General instructions:
    
        [ INITIALIZE ]
        ==============================================================
        python manage.py [opt] [scheme]
        [opt] : Initialize (= Configure + compile)  the jobs.
            opt : 
                --score, --mcore, --acc
                --score : Single core job
                --mcore : Multicore job
                --acc   : GPU accelerated job
            scheme:
                fv : Simulation using finite volume method with 7th order WENO.
                fd : Simulation using finite difference method with 3rd order Kreiss-Oliger dissipation.
	

        eg:
            Run finite volume scheme     : $ python manage.py  --acc fv 
            Run finite difference scheme : $ python manage.py  --acc fd


        [ RUN ]
        ==============================================================

        To run the code after this, traverse to the job directory and run the 
        executable by 
    
        $./main --id <ID> --conf job.config 

        <ID> here can be anythin to tag the output. name of the job directory is preffered.

                            OR

        To restart a truncated job from stored data(in the .bin files), use

        $./main --id <ID> --ff --conf job.config.

        Here <ID> has to be the <ID> used in the truncated job.
	"""
# --------------------------------------------------------------------------------------------------

presets_filename = "presets.hpp"

# --------------------------------------------------------------------------------------------------

TARGET = "main"
OBJECT = "main.o"
SOURCE = "main.cpp"

# --------------------------------------------------------------------------------------------------

GCC = "g++"
PCC = "pgc++"
OPT = "-fast -O3"
STD = "-std=c++0x"

# --------------------------------------------------------------------------------------------------

file_exists_warning_prompt = ("[ WARNING ]...The existing directrories and will be removed. " + 
                            "Do you want to continue? (y/n): "
                              )
# --------------------------------------------------------------------------------------------------

proj_dir = os.getcwd()
resources_dir = os.path.join(proj_dir, "lib")
presets_file = os.path.join(resources_dir, presets_filename)

# Configuration file for the collection of jobs.
batch_configs_file = os.path.join(resources_dir, "configs.yaml")

#List of the jobs will be stored in this txt file.
jobs_list_file = "jobs_list.txt"

# A job.config file will be added to each job folder
job_config_file = ""

# For submitting with condor
condor_conf = os.path.join(resources_dir, "condor.conf")

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
        print(f'{path} removed.')
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
        if config["collective_osc_on"] == True:
            presets.write("#define COLL_OSC_ON\n")
        if config["advection_off"] == True:
            presets.write("#define ADVEC_OFF\n")
        presets.write("\n")
        
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
                N_ITER = int(end_time/dt)+1

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
                    'gz': gz,
                    'v0': v0,
                    'v1': v1,
                    'nz': nz,
                    'N_ITER': N_ITER,
                    'ANAL_EVERY': ANAL_EVERY,
                    'pmo': config['pmo'],
                    'omega': config['omega'],
                    'theta': config['theta'],
                    'mu': config['mu'],
                    'n_fullsnap' : config['n_fullsnap'],
                    'n_vsnap'  : config['n_vsnap'],
                    'vsnap_z'  : config['vsnap_z'],
                    'n_zsnap'  : config['n_zsnap'],
                    'zsnap_v'  : config['zsnap_v'],
                    'n_dump_rho' : config['n_dump_rho'],
                    'dump_rho_v_modes' : config['dump_rho_v_modes'],
                }

                # ID for the job
                config_id = f"{nz}_{nvz}_{CFL}"

                # Making directory with name ID
                model_dir = os.path.join(scheme_dir, config_id)

                # Config file for the job "ID"
                config_file = f"job.config"

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
# --------------------------------------------------------------------------------------------------

def compi(comp_opt = "--acc"):
    """To compile the code"""
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
        comp_stat = os.system(f"make")
        # comp_stat = os.system(f"{PCC} {OPT} -ta=multicore -Minfo=accel -o {TARGET} {SOURCE}")
        if comp_stat == 0:
            print(f"{TARGET} generated for multicore job")
        else:
            return None, None
    elif comp_opt == "--acc":
        comp_stat = os.system(f"{PCC} {OPT} -acc -Minfo=accel -ta=tesla:managed -o {TARGET} {SOURCE}")

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
                
    return "success", jobs_list
    
# --------------------------------------------------------------------------------------------------

def run(jobs_list, scheme_dir_path):
    pwd = os.getcwd()
    for job in jobs_list:
        os.chdir(os.path.join(scheme_dir_path, job))
        os.system(f"./{TARGET} --id {job} --conf job.config")
        os.chdir(pwd)
    return "success"

# --------------------------------------------------------------------------------------------------

def main(mode, scheme):
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
    
##Uncomment this section to run the jobs automatically.
#     if cp_stat == "success":
#         """Run jobs"""
#         run_stat = run(jobs_list, scheme_dir_path)
#     else:
#         print(f"Copying executable failed.")
#         return
        
#     if run_stat == "success":
#         print("SUCCESS")
#     else:
#         print("FAILED")

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    PWD = os.getcwd()
    scheme = None
    mode = None
    is_scheme = False
    is_mode = False
    
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv)):
            if sys.argv[i] == "--help":
                print(instructions)
            elif sys.argv[i] == "--score" or sys.argv[i] == "--mcore" or sys.argv[i] == "--acc":
                mode  = sys.argv[i]                
                i += 1
                is_mode = True
            elif sys.argv[i] == 'fv' or sys.argv[i] == 'fd':
                scheme = sys.argv[i]
                i += 1
                is_scheme = False
            else:
                print(f"Unrecognized option {sys.argv[i]}")
                print(instructions)
                sys.exit()

    if (not is_scheme and not is_mode):
        print(instructions)
    else:
        main(mode, scheme)
        os.chdir(PWD)
