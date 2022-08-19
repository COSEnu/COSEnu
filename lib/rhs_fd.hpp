void NuOsc::calRHS(FieldVar * __restrict out, const FieldVar * __restrict in)
{
    #pragma omp parallel for
    #pragma acc parallel loop
    for (int j = 0; j < nz; j++) {
        
        // common integral over vz'
        real idv_bexR_m_exR  = 0;
        real idv_bexI_p_exI  = 0;
        real ivdv_bexR_m_exR = 0;
        real ivdv_bexI_p_exI = 0;
        real idv_bxx_m_bee_m_xx_p_ee  = 0;
        real ivdv_bxx_m_bee_m_xx_p_ee = 0;

        // OMP reduction not useful here
        #pragma acc loop reduction(+:idv_bexR_m_exR,idv_bexI_p_exI,idv_bxx_m_bee_m_xx_p_ee,ivdv_bexR_m_exR,ivdv_bexI_p_exI,ivdv_bxx_m_bee_m_xx_p_ee)
        for (int k=0;k<nvz; k++) {
            uint kj = idx(k,j);
            idv_bexR_m_exR  +=        (in->bex_re[kj] - in->ex_re[kj] );
            idv_bexI_p_exI  +=        (in->bex_im[kj] + in->ex_im[kj] );
            ivdv_bexR_m_exR +=  vz[k]*(in->bex_re[kj] - in->ex_re[kj] );
            ivdv_bexI_p_exI +=  vz[k]*(in->bex_im[kj] + in->ex_im[kj] );
            idv_bxx_m_bee_m_xx_p_ee  +=       (in->bxx[kj]-in->bee[kj]+in->ee[kj]-in->xx[kj] );
            ivdv_bxx_m_bee_m_xx_p_ee += vz[k]*(in->bxx[kj]-in->bee[kj]+in->ee[kj]-in->xx[kj] );
        }

        // OMP reduction not useful here
        #pragma acc loop vector
        for (int i = 0; i < nvz; i++)
        {
            real *ee = &(in->ee[idx(i, j)]);
            real *xx = &(in->xx[idx(i, j)]);
            real *exr = &(in->ex_re[idx(i, j)]);
            real *exi = &(in->ex_im[idx(i, j)]);
            
            real *bee = &(in->bee[idx(i, j)]);
            real *bxx = &(in->bxx[idx(i, j)]);
            real *bexr = &(in->bex_re[idx(i, j)]);
            real *bexi = &(in->bex_im[idx(i, j)]);

#ifndef ADVEC_OFF
#if defined(ADVEC_CENTER_FD)
            // 2) advection term: (central FD)
            //   4-th order FD for 1st-derivation ~~ ( (a[-2]-a[2])/12 - 2/3*( a[-1]-a[1]) ) / dx
            real factor = -vz[i] / (12 * dz);
#define ADV_FD(x) ((x[-2] - x[2]) - 8.0 * (x[-1] - x[1]))
#elif defined(ADVEC_UPWIND)
            int sv = sgn(vz[i]);
            real factor = -sv * vz[i] / (12 * dz);
#define ADV_FD(x) (-4 * x[-3 * sv] + 18 * x[-2 * sv] - 36 * x[-sv] + 22 * x[0])
#else
            // advection term: (4-th order lopsided finite differencing)
            int sv = sgn(vz[i]);
            real factor = sv * vz[i] / (12 * dz);
#define ADV_FD(x) (x[-3 * sv] - 6 * x[-2 * sv] + 18 * x[-sv] - 10 * x[0] - 3 * x[sv])
#endif
#endif

#ifdef COLL_OSC_ON
            // prepare vz-integral with a cubature rule
            real Iee    = 2*mu* (         exr[0]  *(idv_bexI_p_exI - vz[i]*ivdv_bexI_p_exI ) +  exi[0]*(idv_bexR_m_exR          - vz[i]*ivdv_bexR_m_exR )  );
            real Iexr   =   mu* (   (xx[0]-ee[0]) *(idv_bexI_p_exI - vz[i]*ivdv_bexI_p_exI ) +  exi[0]*(idv_bxx_m_bee_m_xx_p_ee - vz[i]*ivdv_bxx_m_bee_m_xx_p_ee) );
            real Iexi   =   mu* (   (xx[0]-ee[0]) *(idv_bexR_m_exR - vz[i]*ivdv_bexR_m_exR ) -  exr[0]*(idv_bxx_m_bee_m_xx_p_ee - vz[i]*ivdv_bxx_m_bee_m_xx_p_ee) );
            real Ibee   = 2*mu* (        bexr[0]  *(idv_bexI_p_exI - vz[i]*ivdv_bexI_p_exI ) - bexi[0]*(idv_bexR_m_exR          - vz[i]*ivdv_bexR_m_exR )  );
            real Ibexr  =   mu* ( (bxx[0]-bee[0]) *(idv_bexI_p_exI - vz[i]*ivdv_bexI_p_exI ) - bexi[0]*(idv_bxx_m_bee_m_xx_p_ee - vz[i]*ivdv_bxx_m_bee_m_xx_p_ee) );
            real Ibexi  =   mu* ( (bee[0]-bxx[0]) *(idv_bexR_m_exR - vz[i]*ivdv_bexR_m_exR ) + bexr[0]*(idv_bxx_m_bee_m_xx_p_ee - vz[i]*ivdv_bxx_m_bee_m_xx_p_ee) );
#else
            real Iee    = 0;
            real Iexr   = 0;
            real Iexi   = 0;
            real Ibee   = 0;
            real Ibexr  = 0;
            real Ibexi  = 0;
#endif
            // end of mu-part

#ifndef KO_ORD_3
            // Kreiss-Oliger dissipation (5-th order)
            real ko_eps = -ko / dz / 64.0;
#define KO_FD(x) (x[-3] + x[3] - 6 * (x[-2] + x[2]) + 15 * (x[-1] + x[1]) - 20 * x[0])
#else
            // Kreiss-Oliger dissipation (3-nd order)
            real ko_eps = -ko / dz / 16.0;
#define KO_FD(x) (x[-2] + x[2] - 4 * (x[-1] + x[1]) + 6 * x[0])
#endif

#ifdef VAC_OSC_ON
            out->ee[idx(i, j)]     = -pmo * 2 * st * exi[0]                             + ko_eps * KO_FD(ee)   + dv * Iee   + factor * ADV_FD(ee);;
            out->xx[idx(i, j)]     =  pmo * 2 * st * exi[0]                             + ko_eps * KO_FD(xx)   - dv * Iee   + factor * ADV_FD(xx);;
            out->ex_re[idx(i, j)]  = -pmo * 2 * ct * exi[0]                             + ko_eps * KO_FD(exr)  + dv * Iexr  + factor * ADV_FD(exr);;
            out->ex_im[idx(i, j)]  =  pmo * (2 * ct * exr[0] + st * (ee[0] - xx[0]))    + ko_eps * KO_FD(exi)  + dv * Iexi  + factor * ADV_FD(exi);;
            out->bee[idx(i, j)]    = -pmo * 2 * st * bexi[0]                            + ko_eps * KO_FD(bee)  + dv * Ibee  + factor * ADV_FD(bee);;
            out->bxx[idx(i, j)]    =  pmo * 2 * st * bexi[0]                            + ko_eps * KO_FD(bxx)  - dv * Ibee  + factor * ADV_FD(bxx);
            out->bex_re[idx(i, j)] = -pmo * 2 * ct * bexi[0]                            + ko_eps * KO_FD(bexr) + dv * Ibexr + factor * ADV_FD(bexr);;
            out->bex_im[idx(i, j)] =  pmo * (2 * ct * bexr[0] + st * (bee[0] - bxx[0])) + ko_eps * KO_FD(bexi) + dv * Ibexi + factor * ADV_FD(bexi);;
#else
            out->ee[idx(i, j)]     = ko_eps * KO_FD(ee)   + dv * Iee   + factor * ADV_FD(ee);;
            out->xx[idx(i, j)]     = ko_eps * KO_FD(xx)   - dv * Iee   + factor * ADV_FD(xx);;
            out->ex_re[idx(i, j)]  = ko_eps * KO_FD(exr)  + dv * Iexr  + factor * ADV_FD(exr);;
            out->ex_im[idx(i, j)]  = ko_eps * KO_FD(exi)  + dv * Iexi  + factor * ADV_FD(exi);;
            out->bee[idx(i, j)]    = ko_eps * KO_FD(bee)  + dv * Ibee  + factor * ADV_FD(bee);;
            out->bxx[idx(i, j)]    = ko_eps * KO_FD(bxx)  - dv * Ibee  + factor * ADV_FD(bxx);
            out->bex_re[idx(i, j)] = ko_eps * KO_FD(bexr) + dv * Ibexr + factor * ADV_FD(bexr);;
            out->bex_im[idx(i, j)] = ko_eps * KO_FD(bexi) + dv * Ibexi + factor * ADV_FD(bexi);;
#endif

#undef ADV_FD
#undef KO_FD
        }
    }
}

