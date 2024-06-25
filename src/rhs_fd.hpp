void NuOsc::calRHS(FieldVar *out, const FieldVar *in)
{
#pragma omp parallel for collapse(2)
#pragma acc parallel loop  gang /*async*/ // default(present) //gang //private(fac, s, ij, ee, xx, exr, exi, bee, bexx, bexr, bexi)
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

            // 3) interaction terms: vz-integral with a simple trapezoidal rule (can be optimized later)
#ifdef COLL_OSC_ON
            real Iee = 0;
            real Ixx = 0;
            real Iexr = 0;
            real Iexi = 0;
            real Ibee = 0;
            real Ibxx = 0;
            real Ibexr = 0;
            real Ibexi = 0;
            
#pragma acc loop seq /* vector*/ reduction(+:Iee, Ixx, Iexr, Iexi, Ibee, Ibxx, Ibexr, Ibexi) private(eep, xxp, expr, expi, beep, bxxp, bexpr, bexpi) //default(present)
            for (int k = 0; k < nvz; k++)
            { // vz' integral
                real eep = (in->ee[idx(k, j)]);
                real xxp = (in->xx[idx(k, j)]);
                real expr = (in->ex_re[idx(k, j)]);
                real expi = (in->ex_im[idx(k, j)]);
                real beep = (in->bee[idx(k, j)]);
                real bxxp = (in->bxx[idx(k, j)]);
                real bexpr = (in->bex_re[idx(k, j)]);
                real bexpi = (in->bex_im[idx(k, j)]);

                // terms for -i* mu * [rho'-rho_bar', rho]
                Iee += 2 * vw[k] * mu * (1 - vz[i] * vz[k]) * (exr[0] * (expi + bexpi) - exi[0] * (expr - bexpr));
                Ixx += -2 * vw[k] * mu * (1 - vz[i] * vz[k]) * (exr[0] * (expi + bexpi) - exi[0] * (expr - bexpr)); // = -Iee
                Iexr += vw[k] * mu * (1 - vz[i] * vz[k]) * ((xx[0] - ee[0]) * (expi + bexpi) + exi[0] * (eep - xxp - beep + bxxp));
                Iexi += vw[k] * mu * (1 - vz[i] * vz[k]) * (-(xx[0] - ee[0]) * (expr - bexpr) - exr[0] * (eep - xxp - beep + bxxp));
                Ibee += 2 * vw[k] * mu * (1 - vz[i] * vz[k]) * (bexr[0] * (expi + bexpi) + bexi[0] * (expr - bexpr));
                Ibxx += -2 * vw[k] * mu * (1 - vz[i] * vz[k]) * (bexr[0] * (expi + bexpi) + bexi[0] * (expr - bexpr)); // = -Ibee
                Ibexr += vw[k] * mu * (1 - vz[i] * vz[k]) * ((bxx[0] - bee[0]) * (expi + bexpi) - bexi[0] * (eep - xxp - beep + bxxp));
                Ibexi += vw[k] * mu * (1 - vz[i] * vz[k]) * ((bxx[0] - bee[0]) * (expr - bexpr) + bexr[0] * (eep - xxp - beep + bxxp));
            }
            // 3.1) calculate integral with simple trapezoidal rule
            out->ee[idx(i, j)] += dv * Iee;
            out->xx[idx(i, j)] += dv * Ixx;
            out->ex_re[idx(i, j)] += dv * Iexr;
            out->ex_im[idx(i, j)] += dv * Iexi;
            out->bee[idx(i, j)] += dv * Ibee;
            out->bxx[idx(i, j)] += dv * Ibxx;
            out->bex_re[idx(i, j)] += dv * Ibexr;
            out->bex_im[idx(i, j)] += dv * Ibexi;
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

