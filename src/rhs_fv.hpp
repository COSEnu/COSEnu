/*---------------------------------------------------------------------------*/

void NuOsc::WENO7(const double *in, double *rflux, double *lflux)
{
    /*
        7th order WENO reconstruction.
        Inputs:
            u   -> Location about which the reconstruction is performed.
            s   -> s = +1 for forward flux, s = -1 for backward flux.

        Output:
            uul, uur -> Value of the flux at i+1/2 or i-1/2.
    */

    double eps = 1E-6;

    double gamma0r = 4. / 35.;
    double gamma1r = 18. / 35.;
    double gamma2r = 12. / 35.;
    double gamma3r = 1. / 35.;

    double gamma0l = 1. / 35.;
    double gamma1l = 12. / 35.;
    double gamma2l = 18. / 35.;
    double gamma3l = 4. / 35.;

    // Smoothness Indices of the stencils.
    /*
        SI0 -> for stensil S0 = {i, i+1, i+2, i+3} -> r = 0 left shift.
        SI1 -> for stensil S1 = {i-1, i, i+1, i+2} -> r = 1 left shift.
        SI2 -> for stensil S2 = {i-2, i-1, i, i+1} -> r = 2 left shift.
        SI2 -> for stensil S3 = {i-3, i-2, i-1, i} -> r = 3 left shift.
    */
#pragma omp parallel for collapse(2)
#pragma acc parallel loop collapse(2) // default(present)
    // private(u, s, w0r_, w1r_, w2r_, w3r_, w0r, w1r, w2r, w3r, wr, w0l_, w1l_, w2l_, w3l_, w0l, w1l, w2l, w3l, wl,u0r, u1r, u2r, u3r, u0l, u1l, u2l, u3l, SI0, SI1, SI2, SI3)
    for (int i = 0; i < nvz; i++)
    {
        for (int j = -1; j < nz + 1; j++)
        {

            int s = sign(vz[i]);

            const double *u = &in[idx(i, j)];

            double SI0 = u[0] * (2107 * u[0] - 9402 * u[1] + 7042 * u[2] - 1854 * u[3]) + u[1] * (11003 * u[1] - 17246 * u[2] + 4642 * u[3]) + u[2] * (7043 * u[2] - 3882 * u[3]) + 547 * u[3] * u[3];
            double SI1 = u[-1] * (547 * u[-1] - 2522 * u[0] + 1922 * u[1] - 494 * u[2]) + u[0] * (3443 * u[0] - 5966 * u[1] + 1602 * u[2]) + u[1] * (2843 * u[1] - 1642 * u[2]) + 267 * u[2] * u[2];
            double SI2 = u[-2] * (267 * u[-2] - 1642 * u[-1] + 1602 * u[0] - 494 * u[1]) + u[-1] * (2843 * u[-1] - 5966 * u[0] + 1922 * u[1]) + u[0] * (3443 * u[0] - 2522 * u[1]) + 547 * u[1] * u[1];
            double SI3 = u[-3] * (547 * u[-3] - 3882 * u[-2] + 4642 * u[-1] - 1854 * u[0]) + u[-2] * (7043 * u[-2] - 17246 * u[-1] + 7042 * u[0]) + u[-1] * (11003 * u[-1] - 9402 * u[0]) + 2107 * u[0] * u[0];

            // Right going flux
            double w0r_ = gamma0r / pow(eps + SI0, 2);
            double w1r_ = gamma1r / pow(eps + SI1, 2);
            double w2r_ = gamma2r / pow(eps + SI2, 2);
            double w3r_ = gamma3r / pow(eps + SI3, 2);

            double wr = w0r_ + w1r_ + w2r_ + w3r_;

            double w0r = w0r_ / wr;
            double w1r = w1r_ / wr;
            double w2r = w2r_ / wr;
            double w3r = w3r_ / wr;

            double u0r = (1. / 4.) * u[0] + (13. / 12.) * u[1] - (5. / 12.) * u[2] + (1. / 12.) * u[3];       // r = 0
            double u1r = (-1. / 12.) * u[-1] + (7. / 12.) * u[0] + (7. / 12.) * u[1] - (1. / 12.) * u[2];     // r = 1
            double u2r = (1. / 12.) * u[-2] - (5. / 12.) * u[-1] + (13. / 12.) * u[0] + (1. / 4.) * u[1];     // r = 2
            double u3r = (-1. / 4.) * u[-3] + (13. / 12.) * u[-2] - (23. / 12.) * u[-1] + (25. / 12.) * u[0]; // r = 3

            rflux[idx(i, j)] = w0r * u0r + w1r * u1r + w2r * u2r + w3r * u3r;

            // Left going flux
            double w0l_ = gamma0l / pow(eps + SI0, 2);
            double w1l_ = gamma1l / pow(eps + SI1, 2);
            double w2l_ = gamma2l / pow(eps + SI2, 2);
            double w3l_ = gamma3l / pow(eps + SI3, 2);

            double wl = w0l_ + w1l_ + w2l_ + w3l_;

            double w0l = w0l_ / wl;
            double w1l = w1l_ / wl;
            double w2l = w2l_ / wl;
            double w3l = w3l_ / wl;

            double u0l = (-1. / 4.) * u[3] + (13. / 12.) * u[2] - (23. / 12.) * u[1] + (25. / 12.) * u[0]; // r = 0
            double u1l = (1. / 12.) * u[2] - (5. / 12.) * u[1] + (13. / 12.) * u[0] + (1. / 4.) * u[-1];   // r = 1
            double u2l = (-1. / 12.) * u[1] + (7. / 12.) * u[0] + (7. / 12.) * u[-1] - (1. / 12.) * u[-2]; // r = 2
            double u3l = (1. / 4.) * u[0] + (13. / 12.) * u[-1] - (5. / 12.) * u[-2] + (1. / 12.) * u[-3]; // r = 3

            lflux[idx(i, j)] = w0l * u0l + w1l * u1l + w2l * u2l + w3l * u3l;
        }
    }
}
/*---------------------------------------------------------------------------*/

void NuOsc::calFlux(const FieldVar *in, Flux *out)
{
    // right flux
    WENO7(in->ee, out->rflux->ee, out->lflux->ee);
    WENO7(in->xx, out->rflux->xx, out->lflux->xx);
    WENO7(in->ex_re, out->rflux->ex_re, out->lflux->ex_re);
    WENO7(in->ex_im, out->rflux->ex_im, out->lflux->ex_im);

    // Left flux
    WENO7(in->bee, out->rflux->bee, out->lflux->bee);
    WENO7(in->bxx, out->rflux->bxx, out->lflux->bxx);
    WENO7(in->bex_re, out->rflux->bex_re, out->lflux->bex_re);
    WENO7(in->bex_im, out->rflux->bex_im, out->lflux->bex_im);
}

/*---------------------------------------------------------------------------*/

void NuOsc::calRHS(FieldVar *out, const FieldVar *in)
{
    calFlux(in, flux);
    updateBufferZone(flux->rflux);
    updateBufferZone(flux->lflux);

#pragma omp parallel for collapse(2)
#pragma acc parallel loop gang /*async*/ // default(present) //gang //private(fac, s, ij, ee, xx, exr, exi, bee, bexx, bexr, bexi)
    for (int i = 0; i < nvz; i++)
    {
#pragma acc loop vector private(fac, s, ij, ee, xx, exr, exi, bee, bexx, bexr, bexi)
        for (int j = 0; j < nz; j++)
        {
            double *ee = &(in->ee[idx(i, j)]);
            double *xx = &(in->xx[idx(i, j)]);
            double *exr = &(in->ex_re[idx(i, j)]);
            double *exi = &(in->ex_im[idx(i, j)]);

            double *bee = &(in->bee[idx(i, j)]);
            double *bxx = &(in->bxx[idx(i, j)]);
            double *bexr = &(in->bex_re[idx(i, j)]);
            double *bexi = &(in->bex_im[idx(i, j)]);

            out->ee[idx(i, j)] = 0;
            out->xx[idx(i, j)] = 0;
            out->ex_re[idx(i, j)] = 0;
            out->ex_im[idx(i, j)] = 0;
            out->bee[idx(i, j)] = 0;
            out->bxx[idx(i, j)] = 0;
            out->bex_re[idx(i, j)] = 0;
            out->bex_im[idx(i, j)] = 0;

#ifndef ADVEC_OFF

            double fac = -vz[i] / dz;
            int s = sign(vz[i]);
            int ij = idx(i, j);

            out->ee[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->ee[ij] - flux->rflux->ee[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->ee[ij + 1] - flux->lflux->ee[ij]);
            out->xx[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->xx[ij] - flux->rflux->xx[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->xx[ij + 1] - flux->lflux->xx[ij]);
            out->ex_re[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->ex_re[ij] - flux->rflux->ex_re[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->ex_re[ij + 1] - flux->lflux->ex_re[ij]);
            out->ex_im[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->ex_im[ij] - flux->rflux->ex_im[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->ex_im[ij + 1] - flux->lflux->ex_im[ij]);

            out->bee[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->bee[ij] - flux->rflux->bee[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->bee[ij + 1] - flux->lflux->bee[ij]);
            out->bxx[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->bxx[ij] - flux->rflux->bxx[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->bxx[ij + 1] - flux->lflux->bxx[ij]);
            out->bex_re[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->bex_re[ij] - flux->rflux->bex_re[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->bex_re[ij + 1] - flux->lflux->bex_re[ij]);
            out->bex_im[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->bex_im[ij] - flux->rflux->bex_im[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->bex_im[ij + 1] - flux->lflux->bex_im[ij]);
#endif

#ifdef VAC_OSC_ON
            out->ee[idx(i, j)] += -2.0 * st * exi[0] * pmo;
            out->xx[idx(i, j)] += 2.0 * st * exi[0] * pmo;
            out->ex_re[idx(i, j)] += -2.0 * ct * exi[0] * pmo;
            out->ex_im[idx(i, j)] += (2.0 * ct * exr[0] + st * (ee[0] - xx[0])) * pmo;

            out->bee[idx(i, j)] += -0.0;
            out->bxx[idx(i, j)] += 2.0 * st * bexi[0] * pmo;
            out->bex_re[idx(i, j)] += -2.0 * ct * bexi[0] * pmo;
            out->bex_im[idx(i, j)] += (2.0 * ct * bexr[0] + st * (bee[0] - bxx[0])) * pmo;
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

#ifdef COLL_OSC_ON
            double Iee = 0.0;
            double Ixx = 0.0;
            double Iexr = 0.0;
            double Iexi = 0.0;
            double Ibee = 0.0;
            double Ibxx = 0.0;
            double Ibexr = 0.0;
            double Ibexi = 0.0;

            double mut = mu;
#pragma acc loop seq /* vector*/ reduction(+ : Iee, Ixx, Iexr, Iexi, Ibee, Ibxx, Ibexr, Ibexi) private(eep, xxp, expr, expi, beep, bxxp, bexpr, bexpi) // default(present)
            for (int k = 0; k < nvz; k++)
            {

                double eep = (in->ee[idx(k, j)]);
                double xxp = (in->xx[idx(k, j)]);
                double expr = (in->ex_re[idx(k, j)]);
                double expi = (in->ex_im[idx(k, j)]);

                double beep = (in->bee[idx(k, j)]);
                double bxxp = (in->bxx[idx(k, j)]);
                double bexpr = (in->bex_re[idx(k, j)]);
                double bexpi = (in->bex_im[idx(k, j)]);

                Iee += 2.0 * mut * (1.0 - vz[i] * vz[k]) * (exr[0] * (expi + bexpi) - exi[0] * (expr - bexpr));
                Ixx += -2.0 * mut * (1.0 - vz[i] * vz[k]) * (exr[0] * (expi + bexpi) - exi[0] * (expr - bexpr));
                Iexr += 1.0 * mut * (1.0 - vz[i] * vz[k]) * ((xx[0] - ee[0]) * (expi + bexpi) + exi[0] * (eep - xxp - beep + bxxp));
                Iexi += 1.0 * mut * (1.0 - vz[i] * vz[k]) * (-(xx[0] - ee[0]) * (expr - bexpr) - exr[0] * (eep - xxp - beep + bxxp));

                Ibee += 2.0 * mut * (1.0 - vz[i] * vz[k]) * (-bexr[0] * (-expi - bexpi) + bexi[0] * (expr - bexpr));
                Ibxx += -2.0 * mut * (1.0 - vz[i] * vz[k]) * (-bexr[0] * (-expi - bexpi) + bexi[0] * (expr - bexpr));
                Ibexr += 1.0 * mut * (1.0 - vz[i] * vz[k]) * (-(bee[0] - bxx[0]) * (expi + bexpi) - bexi[0] * (eep - xxp - beep + bxxp));
                Ibexi += 1.0 * mut * (1.0 - vz[i] * vz[k]) * ((bee[0] - bxx[0]) * (bexpr - expr) + bexr[0] * (eep - xxp - beep + bxxp));
            }

            out->ee[idx(i, j)] += dv * (Iee);
            out->xx[idx(i, j)] += dv * (Ixx);
            out->ex_re[idx(i, j)] += dv * (Iexr);
            out->ex_im[idx(i, j)] += dv * (Iexi);

            out->bee[idx(i, j)] += dv * (Ibee);
            out->bxx[idx(i, j)] += dv * (Ibxx);
            out->bex_re[idx(i, j)] += dv * (Ibexr);
            out->bex_im[idx(i, j)] += dv * (Ibexi);
#endif
        }
    }
}
