
// +------------------------------------------------------------------------------------------+
// |                                          COSEnu                                          |
// +------------------------------------------------------------------------------------------+
// | Contributors:                                                                            |
// |                                                                                          |
// |	Chun-Yu Lin                                                                           |
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
		else if (strcmp(argv[i], "--ff") == 0)
		{
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
				  << "\t"
				  << "$./main --id <ID> --conf <Configuration file name>\n\n"
				  << "If loading field variables from .bin file:\n"
				  << "\t"
				  << "$./main --id <ID> --conf <Configuration file name> --ff\n"
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

	NuOsc state(pars.z0, pars.z1, pars.nz, pars.nvz, pars.CFL, pars.gz, ID, pars.SCHEME, pars.perturbation_size);
	N_ITER = pars.N_ITER;
	std::cout << std::setw(30) << "NUMBER OFITERATIONS: " << N_ITER << std::endl;

#ifdef VAC_OSC_ON
	// Vacuum oscillation parameters
	state.set_vac_pars(pars.pmo, pars.omega, pars.theta);
#endif
#ifdef COLL_OSC_ON
	// Collective oscillation parameters
	state.set_collective_pars(pars.mu);
#endif

	//......................... INITIALIZING STATE ......................... //

	FieldVar *v_stat0 = new FieldVar(state.size);

	if (!is_ff)
	{
		/*
			Initialize the field variables and angular profiles from the subroutine if is_ff==false.
		*/
		state.initialize();
		state.copy_state(state.v_stat, v_stat0);
		state.write_state0(v_stat0);
		t0 = 1;
	}
	else
	{
		/*
			Initialize the field variables and angular profiles from the binary files if is_ff==true.
		*/
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

	if (!is_ff)
	{
		state.surv_prob(state.v_stat, v_stat0, 0);

		state.dump_rho(state.v_stat, 0);
	}

	// ......................... EVOLVING THE STATE ......................... //

	auto start = std::chrono::steady_clock::now();

	for (int t_iter = t0; t_iter < N_ITER; t_iter++)
	{
		state.step_rk4();

		// ......................... Analysis ......................... //

		if (t_iter % pars.ANAL_EVERY == 0)
		{
#ifdef COLL_OSC_ON
			// Analyze to note the deviation of conserved quantities.
			state.analyse(state.v_stat, P0, 0, t_iter);
#endif
			// Estimates the survival probabilities.
			state.surv_prob(state.v_stat, v_stat0, t_iter);
		}
		// std::cout << "dump_interval: " << pars.dump_interval << "\n";

		if ((t_iter % pars.dump_interval) == 0)
		{
			state.dump_rho(state.v_stat, t_iter);
		}

		if (t_iter % ((int)(N_ITER) / 10) == 0)
		{
			// Write the state of the field variable to a binary file
			// so that the execution of the simulation can be restarted
			// from the last stored values (using the --ff flag)
			// if there is any kind of aboting.
			// state.write_state(t);

			std::cout << " " << std::setprecision(4)
					  << (int)(t_iter * 100.0 / (N_ITER - 1)) << " %"
					  << std::endl;
		}
	}

#ifdef COLL_OSC_ON
	// Estimate the angular distribution of the surviuval probabilities
	// at the end of the simulation.
	delete P0;
#endif

	delete v_stat0;

	std::cout << " 100 %" << std::endl
			  << std::endl
			  << "SIMULATION COMPLETED"
			  << std::endl;

	std::ofstream xtf;
	xtf.open(ID + "time.txt", std::ofstream::out | std::ofstream::trunc);
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	xtf << ID << "\t" << elapsed_seconds.count() << std::endl;
	xtf.close();
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	return (0);
}

// ......................... END OF SIMULATION ......................... //
