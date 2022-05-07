
/*---------------------------------------------------------------------------*/

void NuOsc::write_state(unsigned int t)
{
    std::string file_name = ID + "_state.bin";
    std::string file_name_temp = ID + "_state" + std::to_string(t) + ".bin";
    std::ofstream fpw(file_name_temp, std::ios::out | std::ios::binary);
    long int buffer_size = sizeof(double);

    if (!fpw)
    {
        std::cout << "Failed call to open " << file_name_temp << " in NuOsc::write_state. \n";
        return;
    }

    fpw.write((char *)&t, sizeof(int)); // First entry is the current iteration
    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            int ij = idx(i, j);
            fpw.write((char *)&v_stat->ee[ij], buffer_size);
            fpw.write((char *)&v_stat->xx[ij], buffer_size);
            fpw.write((char *)&v_stat->ex_re[ij], buffer_size);
            fpw.write((char *)&v_stat->ex_im[ij], buffer_size);

            fpw.write((char *)&v_stat->bee[ij], buffer_size);
            fpw.write((char *)&v_stat->bxx[ij], buffer_size);
            fpw.write((char *)&v_stat->bex_re[ij], buffer_size);
            fpw.write((char *)&v_stat->bex_im[ij], buffer_size);
        }
    }
    fpw.close();

    std::rename(file_name_temp.c_str(), file_name.c_str());
}

/*---------------------------------------------------------------------------*/

void NuOsc::write_state0(const FieldVar *stat0)
{
    std::string file_name = ID + "_state0.bin";
    std::ofstream fpw(file_name, std::ios::out | std::ios::binary);
    long int buffer_size = sizeof(double);

    if (!fpw)
    {
        std::cout << "Failed call to open " << file_name << " in NuOsc::write_state0. \n";
    }

    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            int ij = idx(i, j);
            fpw.write((char *)&stat0->ee[ij], buffer_size);
            fpw.write((char *)&stat0->xx[ij], buffer_size);
            fpw.write((char *)&stat0->ex_re[ij], buffer_size);
            fpw.write((char *)&stat0->ex_im[ij], buffer_size);

            fpw.write((char *)&stat0->bee[ij], buffer_size);
            fpw.write((char *)&stat0->bxx[ij], buffer_size);
            fpw.write((char *)&stat0->bex_re[ij], buffer_size);
            fpw.write((char *)&stat0->bex_im[ij], buffer_size);
        }
    }
    fpw.close();
}

/*---------------------------------------------------------------------------*/

int NuOsc::read_state()
{
    std::string file_name = ID + "_state.bin";
    bool is_loading_0 = false;
    if (!file_exists(file_name))
    {
        std::cout << "Unable to find " << file_name << "\n";
        file_name = ID + "_state0.bin";
        is_loading_0 = true;
        std::cout << "Loading from " << file_name << ". \n";
    }

    std::ifstream fpr(file_name, std::ios::in | std::ios::binary);
    if (!fpr)
    {
        std::cout << file_name << " did not open. Exiting.\n";
        exit(1);
    }

    int t;
    if (!is_loading_0)
    {
        fpr.read((char *)&t, sizeof(int));
    }
    else
    {
        t = 1;
    }

    std::cout << "Loading from " << file_name << ". \n";
    long int buffer_size = sizeof(double);
    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            int ij = idx(i, j);
            fpr.read((char *)&v_stat->ee[ij], buffer_size);
            fpr.read((char *)&v_stat->xx[ij], buffer_size);
            fpr.read((char *)&v_stat->ex_re[ij], buffer_size);
            fpr.read((char *)&v_stat->ex_im[ij], buffer_size);

            fpr.read((char *)&v_stat->bee[ij], buffer_size);
            fpr.read((char *)&v_stat->bxx[ij], buffer_size);
            fpr.read((char *)&v_stat->bex_re[ij], buffer_size);
            fpr.read((char *)&v_stat->bex_im[ij], buffer_size);
        }
    }
    fpr.close();
    return t;
}

/*---------------------------------------------------------------------------*/

void NuOsc::read_state0(FieldVar *stat0)
{
    std::string file_name = ID + "_state0.bin";
    std::ifstream fpr(file_name, std::ios::in | std::ios::binary);
    if (!fpr)
    {
        std::cout << file_name << " does not exist. Exiting.\n";
        exit(1);
    }

    long int buffer_size = sizeof(double);
    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            int ij = idx(i, j);
            fpr.read((char *)&stat0->ee[ij], buffer_size);
            fpr.read((char *)&stat0->xx[ij], buffer_size);
            fpr.read((char *)&stat0->ex_re[ij], buffer_size);
            fpr.read((char *)&stat0->ex_im[ij], buffer_size);

            fpr.read((char *)&stat0->bee[ij], buffer_size);
            fpr.read((char *)&stat0->bxx[ij], buffer_size);
            fpr.read((char *)&stat0->bex_re[ij], buffer_size);
            fpr.read((char *)&stat0->bex_im[ij], buffer_size);
        }
    }
    fpr.close();
    return;
}

/*---------------------------------------------------------------------------*/

void NuOsc::read_G0()
{
    std::string file_name = ID + "_G0.bin";
    std::ifstream g_state_in(file_name, std::ios::in | std::ios::binary);
    if (!g_state_in)
    {
        std::cout << file_name << " does not exist. Exiting.\n";
        exit(EXIT_FAILURE);
    }

    long int buffer_size = sizeof(double);
    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            int ij = idx(i, j);
            g_state_in.read((char *)&G0->G[ij], buffer_size);
            g_state_in.read((char *)&G0->bG[ij], buffer_size);
        }
    }
    g_state_in.close();
    return;
}

/*---------------------------------------------------------------------------*/
void NuOsc::output_vsnap(double loc, const int t)
{
    // Sanity check
    if (!(loc >= z0 && loc <= z1))
    {
        std::cout << "Runtime Warning:" << std::endl;
        std::cout << "Inside output_vsnap: "
                  << "Expected location (" << z0 << " <= "
                  << "loc"
                  << " <= " << z1 << ")"
                  << ", given loc = " << loc << std::endl;
        return;
    }

    std::string l = std::to_string(loc);
    l = l.substr(0, l.find_first_of(".") + 3);

    std::string v_snap_fname = ID + "_vsnap_t_" + std::to_string(t) + "_z_" + l + "_.dat";
    std::ofstream v_snap_ofstream;

    v_snap_ofstream.open(v_snap_fname, std::ofstream::out | std::ofstream::trunc);
    if (!v_snap_ofstream)
    {
        std::cout << "Runtime Warning:" << std::endl
                  << "Could not create " << v_snap_fname << "."
                  << std::endl;
        return;
    }
    int j = (int)((loc - z0) / dz - 0.5);
    v_snap_ofstream << "# Phase-space snapshot at (z = " << Z[j] << ", t = " << t << ")"
                    << std::endl;

    for (int i = 0; i < nvz; i++)
    {
        int ij = idx(i, j);
        v_snap_ofstream << std::scientific
                        << vz[i] << "\t"
                        << v_stat->ee[ij] << "\t"
                        << v_stat->xx[ij] << "\t"
                        << v_stat->ex_re[ij] << "\t"
                        << v_stat->ex_im[ij] << "\t"

                        << v_stat->bee[ij] << "\t"
                        << v_stat->bxx[ij] << "\t"
                        << v_stat->bex_re[ij] << "\t"
                        << v_stat->bex_im[ij]
                        << std::endl;
    }
    v_snap_ofstream.close();
}

/*---------------------------------------------------------------------------*/
void NuOsc::output_zsnap(const double vmode, const int t)
{
    // Sanity check
    if (!(vmode >= vz0 && vmode <= vz1))
    {
        std::cout << "Runtime Warning:" << std::endl;
        std::cout << "Inside output_zsnap: "
                  << "Expected v-mode (" << 0 << " <= "
                  << "loc"
                  << " < " << nvz << ")"
                  << ", given v-mode = " << vmode << ": Exiting." << std::endl;
        return;
    }

    int i = (int)((vmode - vz0) / dv - 0.5);
    std::string vstr = std::to_string(vz[i]);

    vstr = vstr.substr(0, vstr.find_first_of(".") + 3);

    std::string z_snap_fname = ID + "_zsnap_t_" + std::to_string(t) + "_v_" + vstr + "_.dat";
    std::ofstream z_snap_ofstream;

    z_snap_ofstream.open(z_snap_fname, std::ofstream::out | std::ofstream::trunc);
    if (!z_snap_ofstream)
    {
        std::cout << "Runtime Error:" << std::endl
                  << "Coul not create " << z_snap_fname << "."
                  << std::endl;
        return;
    }

    z_snap_ofstream << "# Domain snapshot for (vz = " << vz[i] << ", t = " << t << ")" << std::endl;

    for (int j = 0; j < nz; j++)
    {
        int ij = idx(i, j);
        z_snap_ofstream << std::scientific
                        << Z[j] << "\t"
                        << v_stat->ee[ij] << "\t"
                        << v_stat->xx[ij] << "\t"
                        << v_stat->ex_re[ij] << "\t"
                        << v_stat->ex_im[ij] << "\t"

                        << v_stat->bee[ij] << "\t"
                        << v_stat->bxx[ij] << "\t"
                        << v_stat->bex_re[ij] << "\t"
                        << v_stat->bex_im[ij]
                        << std::endl;
    }
    z_snap_ofstream.close();
}
/*---------------------------------------------------------------------------*/

void NuOsc::full_snap(const FieldVar *ivstat, std::string mode)
{
    std::string full_snap_filename = ID + "_v_stat_fullsnap.dat";
    std::ofstream vstat_snap_stream;
    if (mode == "create")
    {
        vstat_snap_stream.open(full_snap_filename, std::ofstream::out | std::ofstream::trunc);
        vstat_snap_stream << "#phys_time = " << phy_time << "\n";
    }
    else if (mode == "app")
    {
        vstat_snap_stream.open(full_snap_filename, std::ofstream::out | std::ofstream::app);
        vstat_snap_stream << "\n" << "#phys_time = " << phy_time << "\n";
    }
    else
    {
        std::cout << "Run time err report from snaps.hpp -> void NuOsc::full_snap\n";
        std::cout << "Unexpected mode: " << mode << " v_stat full snap aborted.\n";
        return;
    }
    if (!vstat_snap_stream)
    {
        std::cout << "Run time err report from snaps.hpp -> void NuOsc::full_snap\n";
        std::cout << "Unable to open " << full_snap_filename << "\n";
        return;
    }
    else
    {
        for (int i = 0; i < nvz; i++)
        {
            for (int j = 0; j < nz; j++)
            {
                int ij = idx(i, j);
                vstat_snap_stream << std::fixed << std::setprecision(20) << std::scientific
                                  << vz[i] << "\t" << Z[j] << "\t"
                                  << 2.0 * ivstat->ex_re[ij] / G0->G[ij] << "\t" << "\t" << -2.0 * ivstat->ex_im[ij] / G0->G[ij] << "\t" << ((ivstat->ee[ij]-ivstat->xx[ij])/ G0->G[ij] ) << "\t"
                                  << 2.0 * ivstat->bex_re[ij] / G0->bG[ij] << "\t" << "\t" << 2.0 * ivstat->bex_im[ij] / G0->bG[ij] << "\t" << ((ivstat->bee[ij]-ivstat->bxx[ij])/ G0->bG[ij] ) << std::endl;
/*dump out rho
                                  << ivstat->ee[ij] << "\t" << ivstat->xx[ij] << "\t" << ivstat->ex_re[ij] << "\t" << ivstat->ex_im[ij] << "\t"
                                  << ivstat->bee[ij] << "\t" << ivstat->bxx[ij] << "\t" << ivstat->bex_re[ij] << "\t" << ivstat->bex_im[ij] << std::endl;
*/
            }
            vstat_snap_stream << std::endl;
        }
        vstat_snap_stream.close();
    }
    return;
}
/*---------------------------------------------------------------------------*/

void NuOsc::dump_rho_v(const FieldVar *ivstate, const double v, std::string fmode)
{
    int vidx = v_idx(v);
    if (!(vidx >= 0 && vidx < nvz))
    {
        std::cout << "Inside NuOsc::dump_rho_v(const FieldVar*, const double , std::string)\n"
                  << "v index = " << vidx << "is out of bound.\n";
        return;
    }
    int prec = (v < 0) ? 5 : 4;

    std::string rdump_filename = ID + "_rho_v_" + std::to_string(v).substr(0, prec) + ".dat";
    std::ofstream rsnap_ofstream;
    if (fmode == "create")
    {
        rsnap_ofstream.open(rdump_filename, std::ofstream::out | std::ofstream::trunc);
    }
    else if (fmode == "app")
    {
        rsnap_ofstream.open(rdump_filename, std::ofstream::out | std::ofstream::app);
    }
    else
    {
        std::cout << "Run time err report from snaps.hpp -> void NuOsc::full_snap\n";
        std::cout << "Unexpected mode: " << fmode << " v_stat full snap aborted.\n";
        return;
    }

    if (!rsnap_ofstream)
    {
        std::cout << "Run time err report from snaps.hpp -> void NuOsc::full_snap\n";
        std::cout << "Unable to open " << rdump_filename << "\n";
        return;
    }
    else
    {
        for (int j = 0; j < nz; j++)
        {
            int ij = idx(vidx, j);
            rsnap_ofstream << Z[j] << "\t"
                           << 2.0 * ivstate->ex_re[ij] / G0->G[ij] << "\t"
                           << -2.0 * ivstate->ex_im[ij] / G0->G[ij] << "\t"
                           << (ivstate->ee[ij] - ivstate->xx[ij]) / G0->G[ij] << "\t"

                           << 2.0 * ivstate->bex_re[ij] / G0->bG[ij] << "\t"
                           << 2.0 * ivstate->bex_im[ij] / G0->bG[ij] << "\t"
                           << (ivstate->bee[ij] - ivstate->bxx[ij]) / G0->bG[ij] << "\n";
        }
        rsnap_ofstream << "\n";
    }
    rsnap_ofstream.close();
}

/*---------------------------------------------------------------------------*/
