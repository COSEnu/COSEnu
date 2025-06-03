#if not defined(__RHS_FD__)
#define __RHS_FD__

/*---------------------------------------------------------------------------*/

void NuOsc::calRHS(FieldVar *out, const FieldVar *in)
{

#ifdef COLL_OSC_ON
            // calculating zeroth moment I and first moment J
    FieldVar *Irho, *Jrho;
    Irho = new FieldVar(nz);
    Jrho = new FieldVar(nz);
#pragma omp parallel 
#pragma acc loop gang
    for(int j=0;j<nz;j++){
        Irho->ee[j]=0.0;
        Irho->xx[j]=0.0;
        Irho->ex_re[j]=0.0;
        Irho->ex_im[j]=0.0;
        Irho->bee[j]=0.0;
        Irho->bxx[j]=0.0;
        Irho->bex_re[j]=0.0;
        Irho->bex_im[j]=0.0;
        Jrho->ee[j]=0.0;
        Jrho->xx[j]=0.0;
        Jrho->ex_re[j]=0.0;
        Jrho->ex_im[j]=0.0;
        Jrho->bee[j]=0.0;
        Jrho->bxx[j]=0.0;
        Jrho->bex_re[j]=0.0;
        Jrho->bex_im[j]=0.0;
        real Iee = 0;
        real Ixx = 0;
        real Iexr = 0;
        real Iexi = 0;
        real Ibee = 0;
        real Ibxx = 0;
        real Ibexr = 0;
        real Ibexi = 0;
        real Jee = 0;
        real Jxx = 0;
        real Jexr = 0;
        real Jexi = 0;
        real Jbee = 0;
        real Jbxx = 0;
        real Jbexr = 0;
        real Jbexi = 0;
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
            Iee += vw[i]*ee[0];
            Ixx += vw[i]*xx[0];
            Iexr += vw[i]*exr[0];
            Iexi += vw[i]*exi[0];
            Ibee += vw[i]*bee[0];
            Ibxx += vw[i]*bxx[0];
            Ibexr += vw[i]*bexr[0];
            Ibexi += vw[i]*bexi[0];
            Jee += vw[i]*vz[i]*ee[0];
            Jxx += vw[i]*vz[i]*xx[0];
            Jexr += vw[i]*vz[i]*exr[0];
            Jexi += vw[i]*vz[i]*exi[0];
            Jbee += vw[i]*vz[i]*bee[0];
            Jbxx += vw[i]*vz[i]*bxx[0];
            Jbexr += vw[i]*vz[i]*bexr[0];
            Jbexi += vw[i]*vz[i]*bexi[0];
        }
        Irho->ee[j]=Iee*dv;
        Irho->xx[j]=Ixx*dv;
        Irho->ex_re[j]=Iexr*dv;
        Irho->ex_im[j]=Iexi*dv;
        Irho->bee[j]=Ibee*dv;
        Irho->bxx[j]=Ibxx*dv;
        Irho->bex_re[j]=Ibexr*dv;
        Irho->bex_im[j]=Ibexi*dv;
        Jrho->ee[j]=Jee*dv;
        Jrho->xx[j]=Jxx*dv;
        Jrho->ex_re[j]=Jexr*dv;
        Jrho->ex_im[j]=Jexi*dv; 
        Jrho->bee[j]=Jbee*dv;
        Jrho->bxx[j]=Jbxx*dv;
        Jrho->bex_re[j]=Jbexr*dv;
        Jrho->bex_im[j]=Jbexi*dv; 
    }
#endif


#pragma omp parallel for collapse(2)
#pragma acc parallel loop gang /*async*/ // default(present) //gang //private(fac, s, ij, ee, xx, exr, exi, bee, bexx, bexr, bexi)
    for (int i = 0; i < nvz; i++)
#pragma acc loop vector private(fac, s, ij, ee, xx, exr, exi, bee, bexx, bexr, bexi)
        for (int j = 0; j < nz; j++)
        {
            real *ee = &(in->ee[idx(i, j)]);
            real *xx = &(in->xx[idx(i, j)]);
            real *exr = &(in->ex_re[idx(i, j)]);
            real *exi = &(in->ex_im[idx(i, j)]);
            real *bee = &(in->bee[idx(i, j)]);
            real *bxx = &(in->bxx[idx(i, j)]);
            real *bexr = &(in->bex_re[idx(i, j)]);
            real *bexi = &(in->bex_im[idx(i, j)]);

            out->ee[idx(i, j)] = 0;
            out->xx[idx(i, j)] = 0;
            out->ex_re[idx(i, j)] = 0;
            out->ex_im[idx(i, j)] = 0;
            out->bee[idx(i, j)] = 0;
            out->bxx[idx(i, j)] = 0;
            out->bex_re[idx(i, j)] = 0;
            out->bex_im[idx(i, j)] = 0;

#ifdef VAC_OSC_ON
            // 1) prepare terms for -i [H0, rho]
            out->ee[idx(i, j)] = -pmo * 2 * st * exi[0];
            out->xx[idx(i, j)] = pmo * 2 * st * exi[0];
            out->ex_re[idx(i, j)] = -pmo * 2 * ct * exi[0];
            out->ex_im[idx(i, j)] = pmo * (2 * ct * exr[0] + st * (ee[0] - xx[0]));
            out->bee[idx(i, j)] = -pmo * 2 * st * bexi[0];
            out->bxx[idx(i, j)] = pmo * 2 * st * bexi[0];
            out->bex_re[idx(i, j)] = -pmo * 2 * ct * bexi[0];
            out->bex_im[idx(i, j)] = pmo * (2 * ct * bexr[0] + st * (bee[0] - bxx[0]));
#endif

#if defined(MAT_OSC_ON)
            out->ee[idx(i, j)] += 0.0;
            out->xx[idx(i, j)] += 0.0;
            out->ex_re[idx(i, j)] += Hm[j] * exi[0];
            out->ex_im[idx(i, j)] += -Hm[j] * exr[0];

            out->bee[idx(i, j)] += 0.0;
            out->bxx[idx(i, j)] += 0.0;
            out->bex_re[idx(i, j)] += -Hm[j] * bexi[0];
            out->bex_im[idx(i, j)] += Hm[j] * bexr[0];
#endif

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
            out->ee[idx(i, j)] += factor * ADV_FD(ee);
            out->xx[idx(i, j)] += factor * ADV_FD(xx);
            out->ex_re[idx(i, j)] += factor * ADV_FD(exr);
            out->ex_im[idx(i, j)] += factor * ADV_FD(exi);
            out->bee[idx(i, j)] += factor * ADV_FD(bee);
            out->bxx[idx(i, j)] += factor * ADV_FD(bxx);
            out->bex_re[idx(i, j)] += factor * ADV_FD(bexr);
            out->bex_im[idx(i, j)] += factor * ADV_FD(bexi);
#undef ADV_FD
#endif

            // 3) interaction terms: 
            // terms for -i* mu * [rho'-rho_bar', rho]

#ifdef COLL_OSC_ON            
            // 3.1) calculate integral with simple trapezoidal rule
            out->ee[idx(i, j)] +=  2*mu* ( (exr[0]*(Irho->ex_im[j]+Irho->bex_im[j])-exi[0]*(Irho->ex_re[j]-Irho->bex_re[j])) 
                                   - vz[i]*(exr[0]*(Jrho->ex_im[j]+Jrho->bex_im[j])-exi[0]*(Jrho->ex_re[j]-Jrho->bex_re[j])) );
            out->xx[idx(i, j)] += -2*mu* ( (exr[0]*(Irho->ex_im[j]+Irho->bex_im[j])-exi[0]*(Irho->ex_re[j]-Irho->bex_re[j])) 
                                   - vz[i]*(exr[0]*(Jrho->ex_im[j]+Jrho->bex_im[j])-exi[0]*(Jrho->ex_re[j]-Jrho->bex_re[j])) );
            out->ex_re[idx(i, j)] += mu* ( ((xx[0]-ee[0])*(Irho->ex_im[j]+Irho->bex_im[j])+exi[0]*(Irho->ee[j]-Irho->xx[j]-Irho->bee[j]+Irho->bxx[j])) 
                                   - vz[i]*((xx[0]-ee[0])*(Jrho->ex_im[j]+Jrho->bex_im[j])+exi[0]*(Jrho->ee[j]-Jrho->xx[j]-Jrho->bee[j]+Jrho->bxx[j])) ); 
            out->ex_im[idx(i, j)] += mu* ( (-(xx[0]-ee[0])*(Irho->ex_re[j]-Irho->bex_re[j])-exr[0]*(Irho->ee[j]-Irho->xx[j]-Irho->bee[j]+Irho->bxx[j])) 
                                   - vz[i]*(-(xx[0]-ee[0])*(Jrho->ex_re[j]-Jrho->bex_re[j])-exr[0]*(Jrho->ee[j]-Jrho->xx[j]-Jrho->bee[j]+Jrho->bxx[j])) ); 
            out->bee[idx(i, j)] +=  2*mu* ( (bexr[0]*(Irho->ex_im[j]+Irho->bex_im[j])+bexi[0]*(Irho->ex_re[j]-Irho->bex_re[j]))  
                                    - vz[i]*(bexr[0]*(Jrho->ex_im[j]+Jrho->bex_im[j])+bexi[0]*(Jrho->ex_re[j]-Jrho->bex_re[j])) ) ;
            out->bxx[idx(i, j)] += -2*mu* ( (bexr[0]*(Irho->ex_im[j]+Irho->bex_im[j])+bexi[0]*(Irho->ex_re[j]-Irho->bex_re[j])) 
                                    - vz[i]*(bexr[0]*(Jrho->ex_im[j]+Jrho->bex_im[j])+bexi[0]*(Jrho->ex_re[j]-Jrho->bex_re[j])) );
            out->bex_re[idx(i, j)] += mu* ( ((bxx[0]-bee[0])*(Irho->ex_im[j]+Irho->bex_im[j])-bexi[0]*(Irho->ee[j]-Irho->xx[j]-Irho->bee[j]+Irho->bxx[j])) 
                                    - vz[i]*((bxx[0]-bee[0])*(Jrho->ex_im[j]+Jrho->bex_im[j])-bexi[0]*(Jrho->ee[j]-Jrho->xx[j]-Jrho->bee[j]+Jrho->bxx[j])) );
            out->bex_im[idx(i, j)] += mu* ( ((bxx[0]-bee[0])*(Irho->ex_re[j]-Irho->bex_re[j])+bexr[0]*(Irho->ee[j]-Irho->xx[j]-Irho->bee[j]+Irho->bxx[j])) 
                                    - vz[i]*((bxx[0]-bee[0])*(Jrho->ex_re[j]-Jrho->bex_re[j])+bexr[0]*(Jrho->ee[j]-Jrho->xx[j]-Jrho->bee[j]+Jrho->bxx[j])) );
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
            out->ee[idx(i, j)] += ko_eps * KO_FD(ee);
            out->xx[idx(i, j)] += ko_eps * KO_FD(xx);
            out->ex_re[idx(i, j)] += ko_eps * KO_FD(exr);
            out->ex_im[idx(i, j)] += ko_eps * KO_FD(exi);
            out->bee[idx(i, j)] += ko_eps * KO_FD(bee);
            out->bxx[idx(i, j)] += ko_eps * KO_FD(bxx);
            out->bex_re[idx(i, j)] += ko_eps * KO_FD(bexr);
            out->bex_im[idx(i, j)] += ko_eps * KO_FD(bexi);
#undef KO_FD
        }
}

/*---------------------------------------------------------------------------*/

#endif // __RHS_FD__

/*---------------------------------------------------------------------------*/
