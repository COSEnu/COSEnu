/*
Initialize the components of \rho and \bar{\rho}
here.
*/
void NuOsc::initialize()
{
    double signu = 0.6;
    double sigbnu = 0.53;
    double alpha = 0.9;
    // Init value
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            G0->G[idx(i, j)] = g(vz[i], 1.0, signu);
            G0->bG[idx(i, j)] = alpha * g(vz[i], 1.0, sigbnu);

            v_stat->ee[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 + eps_(Z[j], 0.0)); 
            v_stat->xx[idx(i, j)]    = 0.5 * G0->G[idx(i, j)] * (1.0 - eps_(Z[j], 0.0));
            v_stat->ex_re[idx(i, j)] = 0.5 * G0->G[idx(i, j)] * (0.0 + eps(Z[j], 0.0));
            v_stat->ex_im[idx(i, j)] = -0.0;

            v_stat->bee[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 + eps_(Z[j], 0.0)); 
            v_stat->bxx[idx(i, j)]    = 0.5 * G0->bG[idx(i, j)] * (1.0 - eps_(Z[j], 0.0));
            v_stat->bex_re[idx(i, j)] = 0.5 * G0->bG[idx(i, j)] * (0.0 + eps(Z[j], 0.0));
            v_stat->bex_im[idx(i, j)] = 0.0;
        }
    }
    updateBufferZone(v_stat);
    std::cout << "Simulation state initialized." << std::endl;
}