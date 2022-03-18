
// +------------------------------------------------------------------------------------------+
// |                                          COSEnu                                          |
// +------------------------------------------------------------------------------------------+
// | Contributors:                                                                            |
// |                                                                                          |
// |	Chun-Yu Lin                                                                            |
// |        - National Center for High-performance computing,                                 |
// |          National Applied Research Laboratories, Hsinchu Science Park,                   |
// |          Hsinchu City 30076, Taiwan.                                                     |
// |                                                                                          |
// |   Meng-Ru Wu, Manu George, Tony Liu, Yi-Siou Wu -                                        |
// |        - Institute of Physics, Academia Sinica, Taipei, 11529, Taiwan.                   |
// |                                                                                          |
// |   Zewei Xiong                                                                            |
// |        - GSI Helmholtzzentrum für Schwerionenforschung, Planckstraße 1, 64291 Darmstadt  |
// |          Germany.                                                                        |
// |                                                                                          |
// |   Kindly direct queries to cosenuproject@gmail.com                                       |
// +------------------------------------------------------------------------------------------+



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
#include <chrono>

using std::cin;
using std::cout;
using std::endl;
using std::string;

#include "presets.hpp"
#include "structures.hpp"
#include "miscellaneous_funcs.hpp"
#include "parser.hpp"
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
	std::string STATE_FILE = "";

	unsigned int N_ITER;
	int t0 = 0;

	bool is_id = false;
	bool is_conf = false;
	bool is_ff = false;

	// ......................... READING COMMANDLINE ARGS ......................... //

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
		else if (strcmp(argv[i], "--ff")==0)
		{
			STATE_FILE = argv[i];
			is_ff = true;
		}
	}

	if (!(is_id && is_conf))
	{
		std::cout << "[ FAIL ]..."
				  << "Both job ID (to label the output files) and configuration should be passed "
				  << "on runtime." << std::endl
				  << "Use --id and --conf to specify them." << std::endl
				  << "General format:\n" 
				  << "\t" << "$./main --id <ID> --conf <Configuration file name>\n"
				  << "If loading field variable from .bin file:\n"
				  << "\t" << "$./main --id <ID> --conf <Configuration file name> --ff\n"
				  << "The --ff flag will look for binary files named ID_state.bin and ID_G0.bin in the same folder.\n"
				  << "eg:(inside condor submit file)" << std::endl
				  << "exiting for now." << std::endl;
		exit(0);
	}
	if (is_ff)
	{
		std::cout << "Initializing from file.\n";
	}

	// ......................... PARSING CONFIG-FILE ......................... //

	Params pars(CONFIG_FILE);

	// ......................... CREATING STATE ......................... //

	NuOsc state(pars.z0, pars.z1, pars.nz, pars.nvz, pars.CFL, pars.gz, ID, pars.SCHEME);
	N_ITER = 1000; //pars.N_ITER;
	std::cout << std::setw(30) << "NUMBER OFITERATIONS: " << N_ITER << std::endl;

#ifdef VAC_OSC_ON
	state.set_vac_pars(pars.pmo, pars.omega, pars.theta);
#endif
#ifdef COLL_OSC_ON
	state.set_collective_pars(pars.mu);
#endif

	//......................... INITIALIZING STATE ......................... //

	FieldVar *v_stat0 = new FieldVar(state.size);

	if (!is_ff)
	{
		state.initialize();
		state.copy_state(state.v_stat, v_stat0);
		state.write_state0(v_stat0);
		t0  = 1;
	}
	else
	{
		state.read_G0();
		t0 = state.read_state();
		state.read_state0(v_stat0);
	}
	std::cout << "Starting time = " << t0 << "\n";


	//......................... EVALUATING INITIAL STATE .........................//

#ifdef COLL_OSC_ON

	Pol *P0 = new Pol(state.size);
	state.cal_pol(v_stat0, P0); // P0 stores initial values of the components of polarization vector.

	if (!is_ff)
	{
		state.analyse(state.v_stat, P0, 0, 0);
	}

#endif


	if(!is_ff)
	{
		state.surv_prob(state.v_stat, v_stat0, 0);
		//............................ SNAPSHOT RLATED ............................//

		for (int i = 0; i < pars.zsnap_v.size(); i++)
		{
			state.output_zsnap(pars.zsnap_v[i], 0);
		}

		// ......................... Phase-space snapshots ......................... //

		for (int i = 0; i < pars.vsnap_z.size(); i++)
		{
			state.output_vsnap(pars.vsnap_z[i], 0);
		}

		// ........................... Full-snapshot ...............................//

		// state.full_snap(state.v_stat, "create");

		// ......................... EVOLVING THE STATE ......................... //

	}

	

    auto start = std::chrono::steady_clock::now();

	for (int t = t0; t < N_ITER; t++)
	{
		state.step_rk4();

		// ...................... Phase-space snapshots ....................... //

		if ((t % pars.v_snap_interval == 0) || (t == N_ITER - 1))
		{
			for (int i = 0; i < pars.vsnap_z.size(); i++)
			{
				state.output_vsnap(pars.vsnap_z[i], t);
			}
		}

		// ......................... Domain snapshots .........................//

		if ((t % pars.z_snap_interval == 0) || (t == N_ITER - 1))
		{
			for (int i = 0; i < pars.zsnap_v.size(); i++)
			{
				state.output_zsnap(pars.zsnap_v[i], t);
			}
		}

		// if ((t % pars.fullsnap_interval) == 0)
		// {
		// 	state.full_snap(state.v_stat, "app");
		// }
		// ......................... Analysis ......................... //
		if (t % pars.ANAL_EVERY == 0)
		{

#ifdef COLL_OSC_ON
			state.analyse(state.v_stat, P0, 0, t);
#endif

			state.surv_prob(state.v_stat, v_stat0, t);
		}

		if (t % ((int)(N_ITER) / 10) == 0)
		{
			state.write_state(t);
			std::cout << " " << std::setprecision(4)
					  << (int)(t * 100.0 / (N_ITER - 1)) << " %"
					  << std::endl;
		}
	}

#ifdef COLL_OSC_ON
	state.v_distr_of_surv_prob(state.v_stat, v_stat0, N_ITER - 1);
	delete P0;
#endif

	delete v_stat0;

	std::cout << " 100 %" << std::endl
			  << std::endl
			  << "SIMULATION COMPLETED"
			  << std::endl;

    std::ofstream xtf;
    xtf.open(ID+"time.txt", std::ofstream::out | std::ofstream::trunc);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    xtf << ID << "\t" <<  elapsed_seconds.count() << std::endl;
    xtf.close();
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    
	return (0);
}

// ......................... END OF SIMULATION ......................... //
