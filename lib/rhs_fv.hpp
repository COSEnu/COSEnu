/*---------------------------------------------------------------------------*/

void NuOsc::WENO7(const double *  __restrict in, double * __restrict rflux, double * __restrict lflux)
{
    /*
        7th order WENO reconstruction.
        Inputs: 
            u   -> Location about which the reconstruction is performed.
            s   -> s = +1 for forward flux, s = -1 for backward flux.

        Output: 
            uul, uur -> Value of the flux at i+1/2 or i-1/2.  
    */

    const double eps = 1E-6;

    const double gamma0r = 4. / 35.;
    const double gamma1r = 18. / 35.;
    const double gamma2r = 12. / 35.;
    const double gamma3r = 1. / 35.;

    const double gamma0l = 1. / 35.;
    const double gamma1l = 12. / 35.;
    const double gamma2l = 18. / 35.;
    const double gamma3l = 4. / 35.;

    const double f1o12 = 1./12.;
    const double f5o12 = 5./12.;
    const double f7o12 = 7./12.;
    const double f13o12 = 13./12.;
    const double f23o12 = 23./12.;
    const double f25o12 = 25./12.;

    //Smoothness Indices of the stencils.
    /*
        SI0 -> for stensil S0 = {i, i+1, i+2, i+3} -> r = 0 left shift.
        SI1 -> for stensil S1 = {i-1, i, i+1, i+2} -> r = 1 left shift.
        SI2 -> for stensil S2 = {i-2, i-1, i, i+1} -> r = 2 left shift.
        SI2 -> for stensil S3 = {i-3, i-2, i-1, i} -> r = 3 left shift.
    */
    #pragma omp parallel for collapse(2)
    #pragma acc parallel loop collapse(2)
    for (int j = -1; j < nz + 1; j++)
    {
        for (int i = 0; i < nvz; i++)
        {

            int s = sign(vz[i]);

            const double *u = &in[idx(i, j)];
            //#pragma acc cache (u[-3*nvz:3*nvz])   // not useful

            double SI0 = u[0     ] * (2107* u[0     ] - 9402 * u[nvz] + 7042 * u[2*nvz] - 1854 * u[3*nvz]) + u[nvz] * (11003 * u[nvz] - 17246 * u[2*nvz] + 4642 * u[3*nvz]) + u[2*nvz] * (7043 * u[2*nvz] - 3882 * u[3*nvz]) + 547 * u[3*nvz] * u[3*nvz];
            double SI1 = u[-  nvz] * (547 * u[-  nvz] - 2522 * u[0] + 1922 * u[nvz] - 494 * u[2*nvz]) + u[0] * (3443 * u[0] - 5966 * u[nvz] + 1602 * u[2*nvz]) + u[nvz] * (2843 * u[nvz] - 1642 * u[2*nvz]) + 267 * u[2*nvz] * u[2*nvz];
            double SI2 = u[-2*nvz] * (267 * u[-2*nvz] - 1642 * u[-nvz] + 1602 * u[0] - 494 * u[nvz]) + u[-nvz] * (2843 * u[-nvz] - 5966 * u[0] + 1922 * u[nvz]) + u[0] * (3443 * u[0] - 2522 * u[nvz]) + 547 * u[1] * u[nvz];
            double SI3 = u[-3*nvz] * (547 * u[-3*nvz] - 3882 * u[-2*nvz] + 4642 * u[-nvz] - 1854 * u[0]) + u[-2*nvz] * (7043 * u[-2*nvz] - 17246 * u[-nvz] + 7042 * u[0]) + u[-1] * (11003 * u[-nvz] - 9402 * u[0]) + 2107 * u[0] * u[0];

            // Right going flux
            double w0r_ = gamma0r / pow(eps + SI0, 2);
            double w1r_ = gamma1r / pow(eps + SI1, 2);
            double w2r_ = gamma2r / pow(eps + SI2, 2);
            double w3r_ = gamma3r / pow(eps + SI3, 2);

            double u0r = ( 0.25  * u[0]    + f13o12 * u[nvz] - f5o12 * u[2*nvz] + f1o12 * u[3*nvz]);       // r = 0
            double u1r = (-f1o12 * u[-nvz] + f7o12  * u[0] + f7o12  * u[nvz] - f1o12 * u[2*nvz]);     // r = 1
            double u2r = ( f1o12 * u[-2*nvz] - f5o12  * u[-nvz] + f13o12 * u[0] + 0.25 * u[nvz]);     // r = 2
            double u3r = (-0.25  * u[-3*nvz] + f13o12 * u[-2*nvz] - f23o12 * u[-nvz] + f25o12 * u[0]); // r = 3

            rflux[idx(i, j)] = ( w0r_ * u0r + w1r_ * u1r + w2r_ * u2r + w3r_ * u3r ) / (w0r_ + w1r_ + w2r_ + w3r_);

            //Left going flux
            double w0l_ = gamma0l / pow(eps + SI0, 2);
            double w1l_ = gamma1l / pow(eps + SI1, 2);
            double w2l_ = gamma2l / pow(eps + SI2, 2);
            double w3l_ = gamma3l / pow(eps + SI3, 2);

            double u0l = -0.25 * u[3*nvz] + f13o12 * u[2*nvz] - f23o12 * u[nvz] + f25o12 * u[0]; // r = 0
            double u1l = f1o12 * u[2*nvz] - f5o12 * u[nvz] + f13o12 * u[0] + 0.25 * u[-nvz];   // r = 1
            double u2l = -f1o12 * u[nvz] + f7o12 * u[0] + f7o12 * u[-nvz] - f1o12* u[-2*nvz]; // r = 2
            double u3l = 0.25 * u[0] + f13o12 * u[-nvz] - f5o12 * u[-2*nvz] + f1o12 * u[-3*nvz]; // r = 3

            lflux[idx(i, j)] = ( w0l_ * u0l + w1l_ * u1l + w2l_ * u2l + w3l_ * u3l ) / (w0l_ + w1l_ + w2l_ + w3l_);
        }
    }
}
/*---------------------------------------------------------------------------*/

void NuOsc::calFlux(const FieldVar * __restrict in, Flux *  __restrict out)
{
    //right flux
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

void NuOsc::calRHS(FieldVar *  __restrict out, const FieldVar * __restrict in)
{
    calFlux(in, flux);
    updateBufferZone(flux->rflux);
    updateBufferZone(flux->lflux);

    #pragma omp parallel for
    #pragma acc parallel loop
    for (int j = 0; j < nz; j++)
    {

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

	// OMP for not useful here
        #pragma acc loop vector
        for (int i = 0; i < nvz; i++)
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

            out->ee[idx(i, j)]    += fac * fabs(1 + s) / 2 * (flux->rflux->ee[ij] - flux->rflux->ee[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->ee[ij + 1] - flux->lflux->ee[ij]);
            out->xx[idx(i, j)]    += fac * fabs(1 + s) / 2 * (flux->rflux->xx[ij] - flux->rflux->xx[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->xx[ij + 1] - flux->lflux->xx[ij]);
            out->ex_re[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->ex_re[ij] - flux->rflux->ex_re[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->ex_re[ij + 1] - flux->lflux->ex_re[ij]);
            out->ex_im[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->ex_im[ij] - flux->rflux->ex_im[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->ex_im[ij + 1] - flux->lflux->ex_im[ij]);
            out->bee[idx(i, j)]    += fac * fabs(1 + s) / 2 * (flux->rflux->bee[ij] - flux->rflux->bee[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->bee[ij + 1] - flux->lflux->bee[ij]);
            out->bxx[idx(i, j)]    += fac * fabs(1 + s) / 2 * (flux->rflux->bxx[ij] - flux->rflux->bxx[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->bxx[ij + 1] - flux->lflux->bxx[ij]);
            out->bex_re[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->bex_re[ij] - flux->rflux->bex_re[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->bex_re[ij + 1] - flux->lflux->bex_re[ij]);
            out->bex_im[idx(i, j)] += fac * fabs(1 + s) / 2 * (flux->rflux->bex_im[ij] - flux->rflux->bex_im[ij - 1]) + fac * fabs(1 - s) / 2 * (flux->lflux->bex_im[ij + 1] - flux->lflux->bex_im[ij]);     
#endif

#ifdef VAC_OSC_ON
            out->ee[idx(i, j)] += -2.0 * st * exi[0] * pmo;
            out->xx[idx(i, j)] += 2.0 * st * exi[0] * pmo;
            out->ex_re[idx(i, j)] += -2.0 * ct * exi[0] * pmo;
            out->ex_im[idx(i, j)] += (2.0 * ct * exr[0] + st * (ee[0] - xx[0])) * pmo;
            out->bee[idx(i, j)] += -2.0 * st * bexi[0] * pmo;
            out->bxx[idx(i, j)] += 2.0 * st * bexi[0] * pmo;
            out->bex_re[idx(i, j)] += -2.0 * ct * bexi[0] * pmo;
            out->bex_im[idx(i, j)] += (2.0 * ct * bexr[0] + st * (bee[0] - bxx[0])) * pmo;
#endif

#ifdef COLL_OSC_ON
            // prepare vz-integral with a cubature rule
            real Iee    = 2*mu* (         exr[0]  *(idv_bexI_p_exI - vz[i]*ivdv_bexI_p_exI ) +  exi[0]*(idv_bexR_m_exR          - vz[i]*ivdv_bexR_m_exR )  );
            real Iexr   =   mu* (   (xx[0]-ee[0]) *(idv_bexI_p_exI - vz[i]*ivdv_bexI_p_exI ) +  exi[0]*(idv_bxx_m_bee_m_xx_p_ee - vz[i]*ivdv_bxx_m_bee_m_xx_p_ee) );
            real Iexi   =   mu* (   (xx[0]-ee[0]) *(idv_bexR_m_exR - vz[i]*ivdv_bexR_m_exR ) -  exr[0]*(idv_bxx_m_bee_m_xx_p_ee - vz[i]*ivdv_bxx_m_bee_m_xx_p_ee) );
            real Ibee   = 2*mu* (        bexr[0]  *(idv_bexI_p_exI - vz[i]*ivdv_bexI_p_exI ) - bexi[0]*(idv_bexR_m_exR          - vz[i]*ivdv_bexR_m_exR )  );
            real Ibexr  =   mu* ( (bxx[0]-bee[0]) *(idv_bexI_p_exI - vz[i]*ivdv_bexI_p_exI ) - bexi[0]*(idv_bxx_m_bee_m_xx_p_ee - vz[i]*ivdv_bxx_m_bee_m_xx_p_ee) );
            real Ibexi  =   mu* ( (bee[0]-bxx[0]) *(idv_bexR_m_exR - vz[i]*ivdv_bexR_m_exR ) + bexr[0]*(idv_bxx_m_bee_m_xx_p_ee - vz[i]*ivdv_bxx_m_bee_m_xx_p_ee) );

            out->ee[idx(i, j)] += dv * (Iee);
            out->xx[idx(i, j)] -= dv * (Iee);
            out->ex_re[idx(i, j)] += dv * (Iexr);
            out->ex_im[idx(i, j)] += dv * (Iexi);
            out->bee[idx(i, j)] += dv * (Ibee);
            out->bxx[idx(i, j)] -= dv * (Ibee);
            out->bex_re[idx(i, j)] += dv * (Ibexr);
            out->bex_im[idx(i, j)] += dv * (Ibexi);
#endif
        }
    }
}

