
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

void NuOsc::dump_rho(const FieldVar *ivstate, const uint t_)
{
    std::string rdump_filename = ID + "_rho_" + std::to_string(t_) + ".dat";
    std::ofstream rsnap_ofstream;
    rsnap_ofstream.open(rdump_filename, std::ofstream::out | std::ofstream::trunc);

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
            for (int i = 0; i < nvz; i++)
            {
                int ij = idx(i, j);
                rsnap_ofstream << Z[j] << "\t" << vz[i] << "\t"
                               << std::setprecision(16)

                               << ivstate->ee[ij] << "\t" << ivstate->xx[ij] << "\t"
                               << ivstate->ex_re[ij] << "\t" << ivstate->ex_im[ij] << "\t"
                               << ivstate->bee[ij] << "\t" << ivstate->bxx[ij] << "\t"
                               << ivstate->bex_re[ij] << "\t" << ivstate->bex_im[ij] << "\n";
            }
            rsnap_ofstream << "\n";
        }
    }
    rsnap_ofstream.close();
}

/*---------------------------------------------------------------------------*/
