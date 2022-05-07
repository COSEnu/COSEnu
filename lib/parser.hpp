/*---------------------------------------------------------------------------*/
template <typename T>
void string_to_type(std::string value, T &var)
{
    /*Convert string to type T and assign it to var*/
    std::stringstream svalue(value);
    svalue >> var;
}

/*---------------------------------------------------------------------------*/

template <typename T>
void cssl_to_vec(std::string cssl, std::vector<T> &vec) // cssl -> comma-separated-string-list
{
    /*
        Convert comma sepparated list of strings(if the fiorm [a, b, c, ...]) and convert it to a 
        vector of type T.
    */
    int start = cssl.find_first_of("[");
    int end = cssl.find_first_of("]");
    cssl = cssl.substr(start + 1, end - 1);

    // In case any accidantal space in the begining of cssl.
    cssl = cssl.substr(cssl.find_first_not_of(" "), end - 1) + ",;";
    std::string c = "";

    while (true)
    {
        c = cssl.substr(0, cssl.find_first_of(","));
        cssl = cssl.substr(cssl.find_first_of(","), cssl.find_first_of("]"));
        cssl = cssl.substr(cssl.find_first_not_of(" ") + 1, cssl.find_first_of("]"));
        std::stringstream svalue(c);
        T tval;
        svalue >> tval;
        vec.push_back(tval);

        if (cssl == ";")
            break;
    }
}

/*---------------------------------------------------------------------------*/

class Params
{
public:
    std::string SCHEME;

    // Simulation specific configs.
    int N_ITER = 0;
    int DUMP_EVERY = 0;
    int ANAL_EVERY = 0;

    // Grid specific configs.
    int nz = 0;
    int nvz = 0;
    int gz = 0;

    double z0 = 0.0;
    double z1 = 0.0;
    double v0;
    double v1;
    double CFL = 0.0;
    double pmo = 0.0;
    double omega = 0.0;
    double theta = 0.0;
    double mu = 1.0;

    //Analysis related
    int n_fullsnap = 0;
    int fullsnap_interval = 0;

    int n_vsnap = 0;
    std::vector<double> vsnap_z;
    int v_snap_interval;

    int n_zsnap = 0;
    std::vector<double> zsnap_v;
    int z_snap_interval;

    int n_dump = 0;
    std::vector<double>v_dumps;
    int dump_interval = 0;


    bool is_scheme = false;
    bool is_nz = false;
    bool is_nvz = false;
    bool is_z0 = false;
    bool is_z1 = false;
    bool is_v0 = false;
    bool is_v1 = false;
    bool is_CFL = false;
    bool is_gz = false;
    bool is_N_ITER = false;
    bool is_ANAL_EVERY = false;
    bool is_pmo = false;
    bool is_omega = false;
    bool is_theta = false;
    bool is_mu = false;

    Params(std::string CONFIG_FILE);
    ~Params() {}
};

Params::Params(std::string CONFIG_FILE)
{
    std::string line;
    std::string key;
    std::string value;

    std::ifstream config;
    config.open(CONFIG_FILE.c_str(), std::ifstream::in);
    if (!config)
    {
        std::cout << "[ FAIL ]...Unable to open "
                  << CONFIG_FILE
                  << "exitting.\n";
        exit(0);
    }
    while (config)
    {
        std::getline(config, line);
        if (line.length() > 0)
        {
            int pos = line.find_first_of(":");
            std::string left = line.substr(0, pos);
            std::string right = line.substr(pos + 1);
            if (left.length() > 0)
            {
                key = left.substr(0, left.find_first_of(" "));
                value = right.substr(right.find_first_not_of(" "), right.find_first_of("\n"));
            }
            if (key == "scheme")
            {
                string_to_type(value, SCHEME);
                is_scheme = true;
            }
            else if (key == "nz")
            {
                string_to_type(value, nz);
                is_nz = true;
            }
            else if (key == "nvz")
            {
                string_to_type(value, nvz);
                is_nvz = true;
            }
            else if (key == "CFL")
            {
                string_to_type(value, CFL);
                is_CFL = true;
            }
            else if (key == "gz")
            {
                string_to_type(value, gz);
                is_gz = true;
            }
            else if (key == "z0")
            {
                string_to_type(value, z0);
                is_z0 = true;
            }
            else if (key == "z1")
            {
                string_to_type(value, z1);
                is_z1 = true;
            }
            else if (key == "v0")
            {
                string_to_type(value, v0);
                is_v0 = true;
            }
            else if (key == "v1")
            {
                string_to_type(value, v1);
                is_v1 = true;
            }
            else if (key == "N_ITER")
            {
                string_to_type(value, N_ITER);
                is_N_ITER = true;
            }
            else if (key == "ANAL_EVERY")
            {
                string_to_type(value, ANAL_EVERY);
                is_ANAL_EVERY = true;
            }
            else if (key == "pmo")
            {
                string_to_type(value, pmo);
                is_pmo = true;
            }
            else if (key == "omega")
            {
                string_to_type(value, omega);
                is_omega = true;
            }
            else if (key == "theta")
            {
                string_to_type(value, theta);
                is_theta = true;
            }
            else if (key == "mu")
            {
                string_to_type(value, mu);
                is_mu = true;
            }
            else if (key == "n_fullsnap")
            {
                string_to_type(value, n_fullsnap);
            }
            else if (key == "n_vsnap")
            {
                string_to_type(value, n_vsnap);
            }
            else if (key == "vsnap_z")
            {
                cssl_to_vec(value, vsnap_z);
            }
            else if (key == "n_zsnap")
            {
                string_to_type(value, n_zsnap);
            }
            else if (key == "zsnap_v")
            {
                cssl_to_vec(value, zsnap_v);
            }
            else if (key == "n_dump_rho")
            {
                string_to_type(value, n_dump);
            }
            else if (key == "dump_rho_v_modes")
            {
                cssl_to_vec(value, v_dumps);
            }

            else
            {
                std::cout << "Redundant or unknown key: " << key << std::endl;
                continue;
            }
        }
    }
    //----------------------------------------------------------

    if (!(is_scheme && is_nz && is_nvz && is_z0 && is_z1 && is_v0 && is_v1 && is_CFL && is_gz &&
          is_N_ITER && is_ANAL_EVERY))
    {
        std::cout << "[ FAIL ]...Incomplete config file" << std::endl
                  << "Exiting." << std::endl;
        exit(0);
    }
    // Snapshot intervals
    // z-snapshot for given vmodes
    if (n_zsnap > 0)
	{
	    z_snap_interval = (int)(N_ITER / n_zsnap);
	}
	else
	{
		z_snap_interval = N_ITER+1;
	}

    // v-snapshots at given locations
    if (n_vsnap > 0)
	{
		v_snap_interval = (int)(N_ITER / n_vsnap);
	}
	else
	{
		v_snap_interval = N_ITER+1;
	}

    // Full snapshots
    if ( n_fullsnap > 0)
	{
	    fullsnap_interval = (int)(N_ITER / n_fullsnap);
            printf("%d %d %d\n",N_ITER,n_fullsnap,fullsnap_interval);
	}
	else
	{
		fullsnap_interval = N_ITER+1;
	}

#ifdef VAC_OSC_ON
    if (!(is_pmo && is_omega && is_theta))
    {
        std::cout << "[ FAIL ]...Incomplete config file" << std::endl
                  << "Exiting." << std::endl;
        exit(0);
    }
#endif
#ifdef COLL_OSC_ON
    if (!(is_mu))
    {
        std::cout << "[ FAIL ]...Incomplete config file" << std::endl
                  << "Exiting." << std::endl;
        exit(0);
    }
    if (n_dump > 0)
    {   
        dump_interval = (int)(N_ITER/n_dump);
    }   
    else
    {   
        dump_interval = N_ITER + 1;
    } 
#endif
}
