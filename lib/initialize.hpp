/*
Initialize the components of \rho and \bar{\rho}
here.
*/
void NuOsc::initialize()
{
    eln_0 = 0;
    beln_0 = 0;
    // Init value
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nvz; i++)
        for (int j = 0; j < nz; j++)
        {
            G0->G[idx(i, j)] = g(vz[i], 1.0, signu);
            G0->bG[idx(i, j)] = alpha * g(vz[i], 1.0, siganu);
        }

// #pragma omp parallel for collapse(2)
    for (int i = 0; i < nvz; i++)
        for (int j = 0; j < nz; j++)
        {
            //step = hev(Z[j], z0/4)-hev(Z[j], z1/4)
            v_stat->ee[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (1.0 + eps_(Z[j], 0.0)); //gauss(Z[j], 0.0, 50);
            v_stat->xx[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (1.0 - eps_(Z[j], 0.0));
            v_stat->ex_re[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + eps(Z[j], 0.0));
            v_stat->ex_im[idx(i, j)] = -0.0;

            eln_0 += (v_stat->ee[idx(i, j)] + v_stat->xx[idx(i, j)])*dz*dv;
        }

// #pragma omp parallel for collapse(2)
    for (int i = 0; i < nvz; i++)
        for (int j = 0; j < nz; j++)
        {
            double th = acos(vz[i]);
            v_stat->bee[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (1.0 + eps_(Z[j], 0.0)); //gauss(Z[j], 0.0, 50);
            v_stat->bxx[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (1.0 - eps_(Z[j], 0.0));
            v_stat->bex_re[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + eps(Z[j], 0.0));
            v_stat->bex_im[idx(i, j)] = 0.0;
            
            beln_0 += (v_stat->bee[idx(i, j)] + v_stat->bxx[idx(i, j)])*dz*dv;
        }
    updateBufferZone(v_stat);
    std::cout << "Simulation state initialized." << std::endl;
}