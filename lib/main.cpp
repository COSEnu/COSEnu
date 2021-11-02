/*
+------------------------------------------------------------------------------------------+
|                                          COSEnu                                          |
+------------------------------------------------------------------------------------------+
| Contributors:                                                                            |
|                                                                                          |
|	Chun-Yu Lin                                                                            |
|        - National Center for High-performance computing,                                 |
|          National Applied Research Laboratories, Hsinchu Science Park,                   |
|          Hsinchu City 30076, Taiwan.                                                     |
|                                                                                          |
|   Meng-Ru Wu, Manu George, Tony Liu, Yi-Siou Wu -                                        |
|        - Institute of Physics, Academia Sinica, Taipei, 11529, Taiwan.                   |
|                                                                                          |
|   Zewei Xiong                                                                            |
|        - GSI Helmholtzzentrum für Schwerionenforschung, Planckstraße 1, 64291 Darmstadt  |
|          Germany.                                                                        |
|                                                                                          |
|   Kindly direct queries to cosenuproject@gmail.com                                       |
+------------------------------------------------------------------------------------------+
*/


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
	bool is_id = false;
	bool is_conf = false;
	unsigned int N_ITER;

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
	}

	if (!(is_id && is_conf))
	{
		std::cout << "[ FAIL ]..."
				  << "Both job ID (to label the output files) and configuration should be passed "
				  << "on runtime." << std::endl
				  << "Use --id and --conf to specify them." << std::endl
				  << "eg:(inside condor submit file)" << std::endl
				  << "arguments = --id $(id) --conf $(conf)" << std::endl
				  << "exiting for now." << std::endl;
		exit(0);
	}

	// ......................... PARSING CONFIG-FILE ......................... //

	Params pars(CONFIG_FILE);

	// ......................... CREATING STATE ......................... //

	NuOsc state(pars.z0, pars.z1, pars.nz, pars.nvz, pars.CFL, pars.gz, ID, pars.SCHEME);
	N_ITER = pars.N_ITER;
	std::cout << std::setw(30) << "NUMBER OFITERATIONS: " << N_ITER << std::endl;

#ifdef VAC_OSC_ON
	state.set_vac_pars(pars.pmo, pars.omega, pars.theta);
#endif
#ifdef COLL_OSC_ON
	state.set_collective_pars(pars.mu);
#endif

	//......................... INITIALIZING STATE ......................... //

	state.initialize();

	//................... MAKING COPYOF THE INITIAL STATE .................. //

	FieldVar *v_stat0 = new FieldVar(state.size);
	state.copy_state(state.v_stat, v_stat0);

	//......................... EVALUATING INITIAL STATE .........................//

#ifdef COLL_OSC_ON
	Pol *P0 = new Pol(state.size);
	state.cal_pol(state.v_stat, P0); // P0 stores initial values of the components of polarization vector.
	state.analyse(state.v_stat, P0, 0, 0);
#endif
	state.survival_prob(state.v_stat, v_stat0, 0);

	//............................ SNAPSHOT RLATED ............................//

	for (int i = 0; i < pars.zsnap_vmodes.size(); i++)
	{
		state.output_zsnap(pars.zsnap_vmodes[i], 0);
	}

	// ......................... Phase-space snapshots ......................... //

	for (int i = 0; i < pars.vsnap_zlocs.size(); i++)
	{
		state.output_vsnap(pars.vsnap_zlocs[i], 0);
	}

	// ........................... Full-snapshot ...............................//

	state.full_snap(state.v_stat, "create");

	// ......................... EVOLVING THE STATE ......................... //

	std::cout << "Running..."
			  << "\n\n";

	for (int t = 1; t < N_ITER; t++)
	{
		state.step_rk4();

		// ...................... Phase-space snapshots ....................... //

		if ((t % pars.v_snap_interval == 0) || (t == N_ITER - 1))
		{
			for (int i = 0; i < pars.vsnap_zlocs.size(); i++)
			{
				state.output_vsnap(pars.vsnap_zlocs[i], t);
			}
		}

		// ......................... Domain snapshots .........................//

		if ((t % pars.z_snap_interval == 0) || (t == N_ITER - 1))
		{
			for (int i = 0; i < pars.zsnap_vmodes.size(); i++)
			{
				state.output_zsnap(pars.zsnap_vmodes[i], t);
			}
		}

		if ((t % pars.fullsnap_interval) == 0)
		{
			state.full_snap(state.v_stat, "app");
		}
		// ......................... Analysis ......................... //
		if (t % pars.ANAL_EVERY == 0)
		{
#ifdef COLL_OSC_ON
			state.analyse(state.v_stat, P0, 0, t);
#endif
			state.survival_prob(state.v_stat, v_stat0, t);
		}

		if (t % ((int)(N_ITER) / 10) == 0)
		{
			std::cout << " " << std::setprecision(4)
					  << (int)(t * 100.0 / (N_ITER - 1)) << " %"
					  << std::endl;
		}
	}

#ifdef COLL_OSC_ON
	state.dom_averaged_survival_prob(state.v_stat, v_stat0, N_ITER - 1);
	delete P0;
#endif
	delete v_stat0;

	std::cout << " 100 %" << std::endl
			  << std::endl
			  << "SIMULATION COMPLETED"
			  << std::endl;

	return (0);
}

// ......................... END OF SIMULATION ......................... //
