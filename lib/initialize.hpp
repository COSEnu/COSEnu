void NuOsc::initialize()
{
/*
    Initialize the components of \rho and \bar{\rho}
    here.
*/
    double signu = 0.6;
    double sigbnu = 0.53;
    double alpha = 0.9;

    int amax=nz/2/10;
    real *phi = new real[nz/10+1];
    real L=z0-z1;
    for(int k=-amax;k<=amax;k++){
        phi[k+amax]=2.0*M_PI*rand()/RAND_MAX;
    }

    std::ofstream g_file(ID + "_G0.bin",std::ofstream::out | std::ofstream::binary);
    if(!g_file)
    {
        std::cout << "Unable to open " << ID+"_G0.bin" << " file from NuOsc::initialise." 
        << "Will not be storing initial angular profiles.\n";
    }

    for (int j=0;j<nz; j++){
       real tmpr=0.0;
       real tmpi=0.0;
       for(int k=-amax;k<amax;k++){
         if(k!=0){
           tmpr=tmpr+1e-7/abs(k)*cos(2.0*M_PI*k*Z[j]/L+phi[k+amax]);
           tmpi=tmpi+1e-7/abs(k)*sin(2.0*M_PI*k*Z[j]/L+phi[k+amax]);
         }
       }
       real tmp2=sqrt(1.0-tmpr*tmpr-tmpi*tmpi);
       for (int i=0;i<nvz; i++){
         G0->G[idx(i, j)] = g(vz[i], 1.0, signu);
         G0->bG[idx(i, j)] = alpha * g(vz[i], 1.0, sigbnu);
         v_stat->ee[idx(i,j)] = 0.5*G0->G[idx(i, j)]*(1.0+tmp2);
         v_stat->xx[idx(i,j)] = 0.5*G0->G[idx(i, j)]*(1.0-tmp2);
         v_stat->ex_re[idx(i,j)] = 0.5*G0->G[idx(i, j)]*tmpr;
         v_stat->ex_im [idx(i,j)] = 0.5*G0->G[idx(i, j)]*tmpi;
         v_stat->bee[idx(i,j)] = 0.5*G0->bG[idx(i, j)]*(1.0+tmp2);
         v_stat->bxx[idx(i,j)] = 0.5*G0->bG[idx(i, j)]*(1.0-tmp2);
         v_stat->bex_re[idx(i,j)]= 0.5*G0->bG[idx(i, j)]*tmpr;
         v_stat->bex_im[idx(i,j)] = -0.5*G0->bG[idx(i, j)]*tmpi;
         g_file.write((char *)&G0->G [idx(i, j)], sizeof(double)); 
         g_file.write((char *)&G0->bG[idx(i, j)], sizeof(double));
       }
    }



/*    for (int i = 0; i < nvz; i++)
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

            g_file.write((char *)&G0->G [idx(i, j)], sizeof(double)); 
            g_file.write((char *)&G0->bG[idx(i, j)], sizeof(double));
        }
    }*/
    updateBufferZone(v_stat);
    std::cout << "Simulation state initialized." << std::endl;
    g_file.close();
}
