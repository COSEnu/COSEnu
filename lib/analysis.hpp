void NuOsc::analyse(const FieldVar *ivstate, const Pol *P0, uint n, uint t)
{
    /****************************************************************************************/
    /* Subroutine to check the deviation of conserved quantities.                           */
    /*                                                                                      */
    /* (1)---                                                                               */
    /*   Estimate maximal deviation of polarization,max{ep_i} for all ep_i belongs to       */
    /*   {P(z_j, v_i) - P0(z_j, v_i)}, i:0->nvz, j:0->nz.                                   */
    /*                                                                                      */
    /* (2)---                                                                               */
    /*   Evaluate M_0                                                                       */
    /****************************************************************************************/
    Pol *P = new Pol(size);
    M *M_0 = new M(0);
    real eELN;
  
    cal_pol(ivstate, P); // Calculating the components of polarization.
    cal_Mn(M_0, P, 0);    // Calculating the components of M_0
    cal_eELN(&eELN,P);
 
    dcon(P, P0, M_0, t, eELN);

    delete P;
    delete M_0;
}

void NuOsc::cal_pol(const FieldVar *inField, Pol *inP)
{
    /* Calculate the polarization vector(\vec {P}) of nu and anu and correspondig |\vec{P}|*/

    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            inP->P1[idx(i, j)] = 2.0 * inField->ex_re[idx(i, j)] / G0->G[idx(i, j)];
            inP->P2[idx(i, j)] = -2.0 * inField->ex_im[idx(i, j)] / G0->G[idx(i, j)];
            inP->P3[idx(i, j)] = (inField->ee[idx(i, j)] - inField->xx[idx(i, j)]) / G0->G[idx(i, j)];
            inP->normP[idx(i, j)] = fabs(sqrt(pow(inP->P1[idx(i, j)], 2) + pow(inP->P2[idx(i, j)], 2) + pow(inP->P3[idx(i, j)], 2)));

            inP->bP1[idx(i, j)] = 2.0 * inField->bex_re[idx(i, j)] / G0->bG[idx(i, j)];
            inP->bP2[idx(i, j)] = 2.0 * inField->bex_im[idx(i, j)] / G0->bG[idx(i, j)];
            inP->bP3[idx(i, j)] = (inField->bee[idx(i, j)] - inField->bxx[idx(i, j)]) / G0->bG[idx(i, j)];
            inP->normbP[idx(i, j)] = fabs(sqrt(pow(inP->bP1[idx(i, j)], 2) + pow(inP->bP2[idx(i, j)], 2) + pow(inP->bP3[idx(i, j)], 2)));
        }
    }
}

/*---------------------------------------------------------------------------*/

void NuOsc::cal_eELN(real *eELN, const Pol *inP)
{
    /*
        Subroutine to calculate the ELN error.
    */
    real nor;
    nor=0.0;
    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            *eELN += 0.5*(G0->G[idx(i, j)] * (1.0-inP->P3[idx(i, j)]) - G0->bG[idx(i, j)] * (1.0-inP->bP3[idx(i, j)])) * dz * dv;
            nor += (G0->G[idx(i, j)] + G0->bG[idx(i, j)])*dz*dv;
        }
    }
    *eELN=*eELN/nor;
}
/*---------------------------------------------------------------------------*/

void NuOsc::cal_Mn(M *inMn, const Pol *inP, unsigned int n)
{
    /*
        Subroutine to calculate the moments. L(vz[i], n) returns the legendre polynomial 
        of order n (up to n=5) for the mode vz[i]
    */
    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            inMn->e1 += L(vz[i], n) * (G0->G[idx(i, j)] * inP->P1[idx(i, j)] - G0->bG[idx(i, j)] * inP->bP1[idx(i, j)]) * dz * dv;
            inMn->e2 += L(vz[i], n) * (G0->G[idx(i, j)] * inP->P2[idx(i, j)] - G0->bG[idx(i, j)] * inP->bP2[idx(i, j)]) * dz * dv;
            inMn->e3 += L(vz[i], n) * (G0->G[idx(i, j)] * inP->P3[idx(i, j)] - G0->bG[idx(i, j)] * inP->bP3[idx(i, j)]) * dz * dv;
        }
    }

    inMn->norm = sqrt(pow(inMn->e1, 2) + pow(inMn->e2, 2) + pow(inMn->e3, 2));
}

/*---------------------------------------------------------------------------*/

void NuOsc::dcon(const Pol *P, const Pol *P0, M *M_0, int t, real eELN)
{
    /*
        Estimates the deviation of conserved quantities from their initial values and write them ot the
        file.
    */
    std::ofstream con_qty_ofstream; // To store deviation of conserved qtys.
    std::string con_qty_fname = ID + "_conserved_quantities.dat";
    if (t == 0)
    {
        con_qty_ofstream.open(con_qty_fname, std::ofstream::out | std::ofstream::trunc);
    }
    else
    {
        con_qty_ofstream.open(con_qty_fname, std::ofstream::out | std::ofstream::app);
    }
    if (!con_qty_ofstream)
    {
        std::cout << "Unable to open " << con_qty_fname << std::endl;
        return;
    }
    else
    {
        if (t == 0)
        {
            con_qty_ofstream << "# [time, dP_max(t), <dP>(t), <dbP(t), |M0|, eELN]\n";
        }
        double dP    = 0.0;
        double dbP   = 0.0;
        double Nee0  = 0.0;
        double Nbee0 = 0.0;
        double avdP  = 0.0;
        double avbdP = 0.0;

        for (int i = 0; i < nvz; i++)
        {
            for (int j = 0; j < nz; j++)
            {
                int ij = idx(i, j);

                dP  = (dP >= fabs(P->normP[ij] - P0->normP[ij])) ? dP : fabs(P->normP[ij] - P0->normP[ij]);
                dbP = (dbP >= fabs(P->normbP[ij] - P0->normbP[ij])) ? dbP : fabs(P->normbP[ij] - P0->normbP[ij]);

                avdP  += fabs(P->normP[ij] - P0->normP[ij]) * G0->G[ij];
                avbdP += fabs(P->normbP[ij] - P0->normbP[ij]) * G0->bG[ij];

                Nee0  += G0->G[ij];
                Nbee0 += G0->bG[ij];
            }
        }
        con_qty_ofstream << phy_time << "\t"
                         << std::fixed << std::setprecision(20)
                         << std::scientific
                         << ((dP >= dbP) ? dP : dbP) << "\t"
                         << avdP / Nee0 << "\t"
                         << avbdP / Nbee0 << "\t"
                         << M_0->norm << "\t" 
                         << eELN << endl;

        con_qty_ofstream.close();
    }
}

/*---------------------------------------------------------------------------*/

void NuOsc::surv_prob(const FieldVar *ivstate, const FieldVar *ivstate0, uint t)
{
    /*
        Total survival probabilities of \nu and \bar\nu.
    */

    std::ofstream surv_prob_ofstream;
    std::string surv_prob_fname = ID + "_survival_probability.dat";
    if (t == 0)
    {
        surv_prob_ofstream.open(surv_prob_fname, std::ofstream::out | std::ofstream::trunc);
        surv_prob_ofstream << "# [time, <Pee>, <Pbee>]" << std::endl;
    }
    else
    {
        surv_prob_ofstream.open(surv_prob_fname, std::ofstream::out | std::ofstream::app);
    }

    if (!surv_prob_ofstream)
    {
        std::cout << "Unable to open " << surv_prob_fname << std::endl;
    }
    else
    {
//        surv_prob_ofstream << "# [time, <Pee>, <Pbee>]" << std::endl;
    }
    double num_Pee   = 0;
    double num_Pbee  = 0;
    double dnom_Pee  = 0;
    double dnom_Pbee = 0;

    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            num_Pee   += (G0->G[idx(i, j)] - ivstate->ee[idx(i, j)]) * dz * dv;
            num_Pbee  += (G0->bG[idx(i, j)] - ivstate->bee[idx(i, j)]) * dz * dv;
            dnom_Pee  += G0->G[idx(i, j)] * dz * dv;
            dnom_Pbee += G0->bG[idx(i, j)] * dz * dv;
//            num_Pee   += ivstate->ee[idx(i, j)] * dz * dv;
//            num_Pbee  += ivstate->bee[idx(i, j)] * dz * dv;
//          dnom_Pee  += G0->G[idx(i, j)] * dz * dv;
//          dnom_Pbee += G0->bG[idx(i, j)] * dz * dv;
//            dnom_Pee  += ivstate0->ee[idx(i, j)] * dz * dv;
//            dnom_Pbee += ivstate0->bee[idx(i, j)] * dz * dv;
        }
    }
    surv_prob_ofstream << phy_time << "\t"
                       << std::setprecision(20)
                       << std::scientific
                       << (dnom_Pee-num_Pee) / dnom_Pee << "\t"
                       << (dnom_Pbee-num_Pbee) / dnom_Pbee << std::endl;
//                       << num_Pee / dnom_Pee << "\t"
//                       << num_Pbee / dnom_Pbee << std::endl;

    surv_prob_ofstream.close();
}

/*---------------------------------------------------------------------------*/

void NuOsc::v_distr_of_surv_prob(const FieldVar *ivstate, const FieldVar *ivstate0, const uint t)
{
    /*
      Survival probabilities of each mode averaged over the domain at a give time.
    */

    std::ofstream av_spv_ofstream;
    std::string av_spv_fname = "dom_avrgd_surv_prob_" + std::to_string(t) + ".dat";
    av_spv_ofstream.open(av_spv_fname, std::ofstream::out | std::ofstream::trunc);
    if (!av_spv_ofstream)
    {
        std::cout << "Unable to open " << av_spv_fname << std::endl;
        return;
    }
    else
    {
        av_spv_ofstream << "# [vz, <P_ee>(v), <P_bee>(v)]" << std::endl;
    }
    for (int i = 0; i < nvz; i++)
    {
        double P_ee   = 0;
        double P_bee  = 0;
        double P_ee0  = 0;
        double P_bee0 = 0;

        for (int j = 0; j < nz; j++)
        {
            P_ee   += ivstate->ee[idx(i, j)] * dz;
            P_bee  += ivstate->bee[idx(i, j)] * dz;
            P_ee0  += ivstate0->ee[idx(i, j)] * dz;
            P_bee0 += ivstate0->bee[idx(i, j)] * dz;
        }
        av_spv_ofstream << std::scientific << vz[i] << "\t"
                        << (P_ee / P_ee0) << "\t"
                        << (P_bee / P_bee0) << std::endl;
    }
    av_spv_ofstream.close();
}

/*------------------------------- EOF ------------------------------------*/

