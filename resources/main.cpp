// ......................... INCLUDES ......................... //

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::string;

#include "presets.hpp"
#include "structures.hpp"
#include "miscellaneous_funcs.hpp"
#include "nuosc.hpp"

#if defined(FD)
#include "rhs_fd.hpp"
#elif defined(FV)
#include "rhs_fv.hpp"
#endif

#include "initialize.hpp"
#include "snaps.hpp"
#include "analysis.hpp"

// ......................... MAIN ......................... //

int main(int argc, char *argv[])
{
	// ......................... VARIABLES ......................... //
	std::string SCHEME = "";
	std::string ID = "";
	std::string CONFIG_FILE = "";
	bool is_id = false;
	bool is_conf = false;
	bool is_scheme = false;

	std::ifstream config;
	std::string line;
	std::string key;
	std::string value;

	// Simulation specific configs.
	int END_TIME = 0;
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
	double sig_nu = 0.0;
	double sig_anu = 0.0;
	double alpha = 0.0;
	double pmo = 0.0;
	double omega = 0.0;
	double theta = 0.0;
	double mu = 1.0;

	//Analysis related
	int n_vsnap = 0;
	std::vector<double> vsnap_zlocs;
	int v_snap_interval;

	int n_zsnap = 0;
	std::vector<int> zsnap_vmodes;
	int z_snap_interval;

	double vmode_P = 0;

	// Some flags.
	bool is_nz = false;
	bool is_nvz = false;
	bool is_z0 = false;
	bool is_z1 = false;
	bool is_v0 = false;
	bool is_v1 = false;
	bool is_CFL = false;
	bool is_gz = false;
	bool is_sig_nu = false;
	bool is_sig_anu = false;
	bool is_alpha = false;
	bool is_END_TIME = false;
	bool is_ANAL_EVERY = false;
	bool is_pmo = false;
	bool is_omega = false;
	bool is_theta = false;
	bool is_mu = false;

	// ......................... READING COMMANDLINE ARGS ......................... //

	// Reading runtime args.
	for (int i = 1; argv[i] != 0; i++)
	{
		if (strcmp(argv[i], "--id") == 0)
		{
			ID = argv[i + 1];
			is_id = true;
			i += 1;
		}

		else if (strcmp(argv[i], "--conf") == 0)
		{
			CONFIG_FILE = argv[i + 1];
			is_conf = true;
			i += 1;
		}
	}

	if (!(is_id && is_conf))
	{
		std::cout << "[ FAIL ]..."
				  << "Both process ID (to label the output files) and configuration should be passed "
				  << "on runtime." << std::endl
				  << "Use --id and --conf to specify them." << std::endl
				  << "eg:(inside condor submit file)" << std::endl
				  << "arguments = --id $(id) --conf $(conf)" << std::endl
				  << "exiting for now." << std::endl;
		exit(0);
	}
	config.open(CONFIG_FILE.c_str(), std::ifstream::in);
	if (!config)
	{
		std::cout << "[ FAIL ]...Unable to open "
				  << CONFIG_FILE
				  << "exitting for now"
				  << std::endl;
		exit(0);
	}

	// ......................... PARSING CONFIG-FILE ......................... //

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
			else if (key == "END_TIME")
			{
				string_to_type(value, END_TIME);
				is_END_TIME = true;
			}
			else if (key == "ANAL_EVERY")
			{
				string_to_type(value, ANAL_EVERY);
				is_ANAL_EVERY = true;
			}
#ifdef VAC_OSC_ON
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
#endif
#ifdef COLL_OSC_ON
			else if (key == "sig_nu")
			{
				string_to_type(value, sig_nu);
				is_sig_nu = true;
			}
			else if (key == "sig_anu")
			{
				string_to_type(value, sig_anu);
				is_sig_anu = true;
			}
			else if (key == "alpha")
			{
				string_to_type(value, alpha);
				is_alpha = true;
			}
			else if (key == "mu")
			{
				string_to_type(value, mu);
				is_mu = true;
			}
			else if (key == "vmode_P")
			{
				string_to_type(value, vmode_P);
			}
#endif
			else if (key == "n_vsnap")
			{
				string_to_type(value, n_vsnap);
			}
			else if (key == "vsnap_zlocs")
			{
				cssl_to_vec(value, vsnap_zlocs);
			}
			else if (key == "n_zsnaps")
			{
				string_to_type(value, n_zsnap);
			}
			else if (key == "zsnap_vmodes")
			{
				cssl_to_vec(value, zsnap_vmodes);
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
		  is_END_TIME && is_ANAL_EVERY))
	{
		std::cout << "[ FAIL ]...Incomplete config file" << std::endl
				  << "Exiting." << std::endl;
		exit(0);
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
	if (!(is_sig_nu && is_sig_anu && is_mu && is_alpha))
	{
		std::cout << "[ FAIL ]...Incomplete config file" << std::endl
				  << "Exiting." << std::endl;
		exit(0);
	}
#endif

	// ......................... CREATING STATE ......................... //

	NuOsc state(z0, z1, nz, nvz, CFL, gz, ID, SCHEME);
	std::cout << std::setw(30) << "END_TIME: " << END_TIME << std::endl;

#ifdef VAC_OSC_ON
	state.set_vac_pars(pmo, omega, theta);
#endif
#ifdef COLL_OSC_ON
	state.set_collective_pars(mu, sig_nu, sig_anu, alpha);
#endif

	//......................... INITIALIZING STATE ......................... //

	state.initialize();
	FieldVar * v_stat0 = new FieldVar(state.size);
	state.copy_state(state.v_stat, v_stat0);
	//......................... OUTPUT FILE STREAMS .........................//

#ifdef COLL_OSC_ON
	std::ofstream con_qty_ofstream; // To store deviation of conserved qtys.
	std::string conserved_fname = ID + "_conserved_quantities.dat";

	std::ofstream av_surv_prob_ofstream;
	std::string av_surv_prob_fname = ID + "_averaged_survival_prob.dat";

	av_surv_prob_ofstream.open(av_surv_prob_fname, std::ofstream::out | std::ofstream::trunc);
	if (!av_surv_prob_ofstream)
	{
		std::cout << "Unable to open " << av_surv_prob_fname << std::endl;
	}
	else
	{
		av_surv_prob_ofstream << "# [time, <Pee>, <Pbee>]" << std::endl;
	}

	con_qty_ofstream.open(conserved_fname, std::ofstream::out | std::ofstream::trunc);
	if (!con_qty_ofstream)
	{
		std::cout << "Unable to open " << conserved_fname << std::endl;
	}
	else
	{
		con_qty_ofstream << "# [time, dP_max(t), <dP>(t), <dbP(t), <dP>(t, v="
						 << std::to_string(vmode_P) << "), |M0|]" << std::endl;
	}
#endif

	std::ofstream surv_prob_ofstream;
	std::string surv_prob_fname = ID + ".total_survival_prob.dat";
	surv_prob_ofstream.open(surv_prob_fname, std::ofstream::out | std::ofstream::trunc);
	if (!surv_prob_ofstream)
	{
		cout << "Unable to open " << surv_prob_fname << endl;
	}

	//......................... EVALUATING INITIAL STATE .........................//

#ifdef COLL_OSC_ON
	Pol *P0 = new Pol(state.size);
	state.cal_P(state.v_stat, P0); // P0 stores initial polarization.
	state.analyse(state.v_stat, P0, con_qty_ofstream, 0, 0);
	state.averaged_survival_prob(state.v_stat, state.G0, av_surv_prob_ofstream, 0);
#endif
	state.survival_prob(state.v_stat, surv_prob_ofstream, 0);

	//......................... INITIAL STATE SNAPSHOTS .........................//

	// ......................... Domain snapshots .........................//
	int t_zsnap = 0;
	if (n_zsnap > 0)
	{
		z_snap_interval = (int)(END_TIME / n_zsnap);
	}
	else
	{
		z_snap_interval = 0;
	}
	for (int i = 0; i < zsnap_vmodes.size(); i++)
	{
		state.output_zsnap(zsnap_vmodes[i], t_zsnap);
	}
	t_zsnap += z_snap_interval;

	// ......................... Phase-space snapshots ......................... //

	int t_vsnap = 0;
	if (n_vsnap > 0)
	{
		v_snap_interval = (int)(END_TIME / n_vsnap);
	}
	else
	{
		v_snap_interval = 0;
	}

	for (int i = 0; i < vsnap_zlocs.size(); i++)
	{
		state.output_vsnap(vsnap_zlocs[i], t_vsnap);
	}
	t_vsnap += v_snap_interval;

#ifdef COLL_OSC_ON
	state.averaged_survival_prob_v(state.v_stat, v_stat0, 0);
#endif
	// ......................... EVOLVING THE STATE ......................... //

	std::cout << "Running..." << std::endl
			  << std::endl;
	for (int t = 1; t < END_TIME; t++)
	{
		state.step_rk4();

		// ......................... Phase-space snapshots ......................... //

		if (t == t_vsnap)
		{
			for (int i = 0; i < vsnap_zlocs.size(); i++)
			{
				state.output_vsnap(vsnap_zlocs[i], t);
			}
			t_vsnap += v_snap_interval;
		}

		// ......................... Domain snapshots .........................//

		if (t == t_zsnap)
		{
			for (int i = 0; i < zsnap_vmodes.size(); i++)
			{
				state.output_zsnap(zsnap_vmodes[i], t);
			}
			t_zsnap += z_snap_interval;
		}

		// ......................... Analysis ......................... //
		if (t % ANAL_EVERY == 0)
		{
#ifdef COLL_OSC_ON
			state.analyse(state.v_stat, P0, con_qty_ofstream, 0, t);
			state.averaged_survival_prob(state.v_stat, state.G0, av_surv_prob_ofstream, t);
#endif
			state.survival_prob(state.v_stat, surv_prob_ofstream, t);
		}
		if (t % 100 == 0)
		{
			std::cout << " " << std::setprecision(4)
					  << ((float)t * 100 / (END_TIME - 1)) << " %"
					  << std::endl;
		}
	}
	std::cout << std::endl;

#ifdef COLL_OSC_ON
	state.averaged_survival_prob_v(state.v_stat, v_stat0, END_TIME-1);
#endif

#ifdef COLL_OSC_ON
	con_qty_ofstream.close();
	av_surv_prob_ofstream.close();
#endif
	surv_prob_ofstream.close();

	std::cout << "SIMULATION COMPLETED"
			  << std::endl;

	delete v_stat0;
	return (0);
}

// ......................... END OF SIMULATION ......................... //
