/*---------------------------------------------------------------------------*/
void NuOsc::output_vsnap(double loc, const int t)
{
    //Sanity check
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
    //Sanity check
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
                  << "Could not create " << z_snap_fname << "."
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
