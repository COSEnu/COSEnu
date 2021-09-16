# --------------------------------------------------------------------------------------------------

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
file_exists_warning_prompt = ("[ WARNING ]...The existing directrories and will be removed. " +
                              "Do you want to continue? (y/n): "
                              )
# --------------------------------------------------------------------------------------------------
"""
Files to be copied from resources folder to the
job folder.
"""
presets_filename = "presets.hpp"
source_files = [
    'Makefile',
    'Makefile.inc',
    'structures.hpp',
    presets_filename,
    'main.cpp',
    'nuosc.hpp',
    'rhs_fv.hpp',
    'rhs_fd.hpp',
    'initialize.hpp',
    'snaps.hpp',
    'analysis.hpp',
    'miscellaneous_funcs.hpp',
]
# --------------------------------------------------------------------------------------------------

TARGET = "main"
OBJECT = "main.o"

# --------------------------------------------------------------------------------------------------

proj_dir = os.getcwd()
resources_dir = os.path.join(proj_dir, "lib")

# Configuration file for the collection of jobs.
batch_configs_file = os.path.join(resources_dir, "configs.yaml")

#List of the jobs will be stored in this txt file.
jobs_list_file = "jobs_list.txt"

# A job.config file will be added to each job folder
job_config_file = ""

# For submitting with condor
job_submit_src_file = os.path.join(resources_dir, "job.submit")

presets_file = os.path.join(resources_dir, presets_filename)

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


def configure(scheme="fv"):

    proj_dir = os.getcwd()

    with open(batch_configs_file) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    nzs = config["nzs"]
    nvzs = config["nvzs"]
    CFLS = config["CFLS"]
    z = config["zrange"]
    v0 = config["v0"]
    v1 = config["v1"]
    sig_nu = config["signu"]
    sig_anu = config["siganu"]
    alpha = config["alpha"]
    ncycle = config["ncycle"]
    nanalyze = config["nanalyze"]
    scheme_dir = os.path.join(proj_dir,  config[f'folder_{scheme}'])

    with open(presets_file, "w") as presets:
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

    # Creating necessary folders and copying files
    config_list = []
    copied_files = []

    continue_if_exist = "N"

    try:
        os.mkdir(scheme_dir)
        print(f"[ OK ]...Creating {scheme_dir}")
    except FileExistsError:
        print(f"[ EXISTS ]...{scheme_dir}")
        continue_if_exist = input(file_exists_warning_prompt).upper()
        if continue_if_exist == "N":
            print(f"[ OK ]...exiting the {scheme} setup.")
            sys.exit()
        elif continue_if_exist == "Y":
            print(f"[ OK ]...Setting up the projects for {scheme}.")
        else:
            print(f"[ FAIL ]...Invalid input {continue_if_exist} Exiting.")
            sys.exit()

    os.chdir(scheme_dir)
    # copying files from resources directory to the scheme directory
    for file in source_files:
        src = os.path.join(resources_dir, file)
        dst = os.path.join(scheme_dir, file)

        if os.path.exists(scheme_dir):
            try:
                shutil.copy(src, dst)
                print(f"[ OK ]...Copying {src} to {dst}")
                copied_files.append(dst)
            except:
                print(f"[ FAIL ]...Copying {src} to {dst} exiting.")
                sys.exit()
        else:
            print(f"[ FAIL ]...Unable to locate {scheme_dir} Exiting.")
            sys.exit()

    if (scheme == 'fd'):
        gz = 2
    elif (scheme == 'fv'):
        gz = 4
    # Loding config values from config file
    for nz in nzs:
        for nvz in nvzs:
            for CFL in CFLS:
                z0 = z[0]
                z1 = z[1]
                dz = (z1-z0)/nz
                dvz = (v0-v1)/nvz
                dt = CFL*dz/v1
                END_TIME = int(ncycle*z1/dt)
                if nanalyze != 0:
                    ANAL_EVERY = int(END_TIME/nanalyze)
                    if (ANAL_EVERY < 1):
                        ANAL_EVERY = 1
                else:
                    ANAL_EVERY = END_TIME + 1

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
                    'dz': dz,
                    'dvz': dvz,
                    'dt': dt,
                    'END_TIME': END_TIME,
                    'ANAL_EVERY': ANAL_EVERY,
                    'sig_nu': sig_nu,
                    'sig_anu': sig_anu,
                    'alpha': alpha,
                    'pmo': config['pmo'],
                    'omega': config['omega'],
                    'theta': config['theta'],
                    'mu': config['mu'],
                    'n_vsnap' : config['n_vsnap'],
                    'vsnap_zlocs' : config['vsnap_zlocs'],
                    'n_zsnaps' : config['n_zsnaps'],
                    'zsnap_vmodes' : config['zsnap_vmodes'],
                    'vmode_P' : config['vmode_P'],
                }

                # ID for the job
                config_id = f"{nz}_{nvz}_{CFL}"

                # Making directory with name ID
                model_dir = config_id

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
    shutil.copy(job_submit_src_file, scheme_dir)

    os.chdir(proj_dir)

    return proj_dir, scheme_dir, job_config_file, copied_files

# --------------------------------------------------------------------------------------------------


def initialize(scheme):

    proj_dir, scheme_dir, job_config_file, copied_files = configure(
        scheme=scheme)
    os.chdir(scheme_dir)

    # Compiling
    TARGET = "main"
    OBJECT = "main.o"

    for f in os.listdir("./"):
        if f == TARGET or f == OBJECT:
            os.remove(f)

    try:
        print()
        stat_make = os.system(f"make")
        print()
        if stat_make == 0:
            print(f"[ OK ]...Successfully compiled.")
        else:
            print(f"[ FAIL ]...Compilation failed.")
            sys.exit()
    except:
        print(f"[ FAIL ]...Compile")

    for file_item in copied_files:
        os.remove(file_item)

    os.chdir(proj_dir)

    return proj_dir, scheme_dir, job_config_file

# --------------------------------------------------------------------------------------------------


def submit_with_condor(scheme):
    print(
        f'[ START ]...Identifying the scheme directory from {batch_configs_file}.')

    with open(batch_configs_file) as f:
        require = yaml.load(f, Loader=yaml.FullLoader)
        scheme_dir = os.path.join(proj_dir, require[f'folder_{scheme}'])

    print(f'[ FINISHED ]...Scheme directory identified to be {scheme_dir}')
    print(f'[ OK ]...Continuing with the submission.')

    try:
        os.chdir(scheme_dir)
        job_config_file = os.path.join(scheme_dir, jobs_list_file)
        submit_jobs_list = []
        with open(job_config_file, 'r') as f:
            for line in f.readlines():
                submit_jobs_list.append(line.strip("\n").split(","))
    except FileNotFoundError:
        print(f'[ FAIL ]...unable to locate {scheme_dir}. Exiting')
        sys.exit()

    for item in submit_jobs_list:
        id, config_file = item[0], item[1]
        model_dir = os.path.join(scheme_dir, id)

        try:
            pwd = os.getcwd()
            os.chdir(model_dir)
            shutil.copy(job_submit_src_file, os.path.join(
                model_dir, 'job.submit'))
            shutil.copy(os.path.join(scheme_dir, TARGET),
                        os.path.join(model_dir, TARGET))
            model_job_config_file = os.path.join(model_dir, f'job.req')

            with open(model_job_config_file, 'w') as f:
                f.write(f"{id}" + ", " + f"{config_file}")

            os.system(f"condor_submit job.submit")
            os.chdir(pwd)

        except FileNotFoundError:
            print(f'[ FAIL ]...unable to locate {model_dir}. Exiting')
            sys.exit()

    os.chdir(scheme_dir)

# --------------------------------------------------------------------------------------------------


def submit_to_local_machine(scheme):
    import time
    print(
        f'[ START ]...Identifying the job directory from {batch_configs_file}.')

    with open(batch_configs_file) as f:
        require = yaml.load(f, Loader=yaml.FullLoader)
        scheme_dir = os.path.join(proj_dir, require[f'folder_{scheme}'])

    print(f'[ FINISHED ]...job directory identified to be {scheme_dir}')
    print(f'[ OK ]...Continuing with the submission.')

    if os.path.exists(scheme_dir):
        os.chdir(scheme_dir)
        job_config_file = os.path.join(scheme_dir, jobs_list_file)
        submit_jobs_list = []
        try:
            with open(job_config_file, 'r') as f:
                for line in f.readlines():
                    submit_jobs_list.append(line.strip("\n").split(","))
        except FileNotFoundError:
            print(f'[ FAIL ]...unable to find {job_config_file}. Exiting')
            sys.exit()
    else:
        print(f'[ FAIL ]...Path {scheme_dir} does not exists.')
        sys.exit()

    for item in submit_jobs_list:
        id, config_file = item[0], item[1]
        model_dir = os.path.join(scheme_dir, id)

        try:
            pwd = os.getcwd()
            os.chdir(model_dir)
            shutil.copy(os.path.join(scheme_dir, TARGET),
                        os.path.join(model_dir, TARGET))

            print(f"Running {id}")
            start = time.time()
            os.system(f"./{TARGET} --id {id} --conf {config_file}")
            end = time.time()
            print(f"Job: {id} finished in {end-start:.5f} seconds.")
            os.chdir(pwd)

        except FileNotFoundError:
            print(f'[ FAIL ]...unable to locate {model_dir}. Exiting')
            sys.exit()

    os.chdir(scheme_dir)
    try:
        os.remove(OBJECT)
    except FileNotFoundError:
        print(f"[ FAIL ]...Unable to remove {OBJECT} as it does not exist.")

# --------------------------------------------------------------------------------------------------


if __name__ == "__main__":

    available_options = """
	[ --in ] : (Initialize = Configure and compile)  the jobs.
	[ --ls ] : Run the initialize jobs on local machine.
	[ --cs ] : Submit the initialized jobs to the cluster.
	"""

    sample_templates = """
	To initialize finite volume scheme     : $ python manage.py --in fv 
	To initialize finite difference scheme : $ python manage.py --in fd
	To run the initialized jobs locally    : $ python manage.py --ls fd(or fv)
	To submit to the cluster using condor  : $ python manage.py --cs fd(or fv) 
	"""

    PWD = os.getcwd()
    init_scheme = None  # Scheme to be compiled and configured
    cluster_sub_scheme = None
    local_run_scheme = None

    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv)):
            if sys.argv[i] == "--in":
                init_scheme = sys.argv[i+1]
                i += 1
            elif sys.argv[i] == '--cs':
                cluster_sub_scheme = sys.argv[i+1]
                i += 1
            elif sys.argv[i] == '--ls':
                local_run_scheme = sys.argv[i+1]
                i += 1

    if init_scheme == None and cluster_sub_scheme == None and local_run_scheme == None:
        print(f'[ JOBLESS ]...Asked to do nothing. Exiting.')
        print(f"Suggestions:")
        print(f"you have the following options: \n {available_options}")
        print(f"Try using from the following templates: \n {sample_templates}")
        os.chdir(PWD)
        sys.exit()

    """
	Initializing the schemes: 
	Making the necessary directories, 
	copying files from resources and 
	compiling the files.
	"""
    if init_scheme != None:
        if init_scheme == "fv" or init_scheme == "fd":
            proj_dir, scheme_dir, job_config_file = initialize(init_scheme)
        else:
            print(f"[ FAIL ]...Unrecognized input after --init {init_scheme}.")
            print(f"Suggestions:")
            print(f"you have the following options: \n {available_options}")
            print(
                f"Try using from the following templates: \n {sample_templates}")
            os.chdir(PWD)
            sys.exit()

    """
	Submiting the initialized jobs to the cluster.
	"""
    if cluster_sub_scheme != None:
        if cluster_sub_scheme == "fv" or cluster_sub_scheme == "fd":
            pwd = os.getcwd()
            submit_with_condor(cluster_sub_scheme)
            os.chdir(pwd)
        else:
            print(
                f"[ FAIL ]...Unrecognized input {cluster_sub_scheme} after --cs")
            print(
                f"Try using from the following templates: \n {sample_templates}")
            os.chdir(PWD)
            sys.exit()

    """
	Submiting the initialized jobs to the local machine.
	"""
    if local_run_scheme != None:
        if local_run_scheme == "fv" or local_run_scheme == "fd":
            pwd = os.getcwd()
            submit_to_local_machine(local_run_scheme)
            os.chdir(pwd)
        else:
            print(
                f"[ FAIL ]...Unrecognized input {local_run_scheme} after --ls")
            print(
                f"Try using from the following templates: \n {sample_templates}")
            os.chdir(PWD)
            sys.exit()

    os.chdir(PWD)
    print()
    print('[ FINISHED ]')

# --------------------------------------------------------------------------------------------------
