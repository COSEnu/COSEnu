void NuOsc::analyse(const FieldVar *inField, const Pol *inP0, std::ofstream &con_qty_ofstream,
                    unsigned int n, const int time)
{

    Pol *P = new Pol(size);
    M *M0 = new M(0);

    cal_P(inField, P); // Calculating the components of polarization.
    cal_Mn(M0, P, 0);  // Calculating the components of M_0

    double dP = 0.0;
    double dbP = 0.0;
    double Nee0 = 0.0;
    double Nbee0 = 0.0;
    double meandP = 0.0;
    double meanbdP = 0.0;
    double avP_v = 0.0;

    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            dP = (dP >= fabs(P->normP[idx(i, j)] - inP0->normP[idx(i, j)])) ? dP : fabs(P->normP[idx(i, j)] - inP0->normP[idx(i, j)]);
            dbP = (dbP >= fabs(P->normbP[idx(i, j)] - inP0->normbP[idx(i, j)])) ? dbP : fabs(P->normbP[idx(i, j)] - inP0->normbP[idx(i, j)]);

            meandP += fabs(P->normP[idx(i, j)] - inP0->normP[idx(i, j)]) * G0->G[idx(i, j)] * dz * dv;
            meanbdP += fabs(P->normbP[idx(i, j)] - inP0->normbP[idx(i, j)]) * G0->bG[idx(i, j)] * dz * dv;
            Nee0 += G0->G[idx(i, j)] * dz * dv;
            Nbee0 += G0->bG[idx(i, j)] * dz * dv;
        }
    }

    con_qty_ofstream << time << "\t"
                     << std::scientific
                     << ((dP >= dbP) ? dP : dbP) << "\t"
                     << meandP / Nee0 << "\t"
                     << meanbdP / Nbee0 << "\t"
                     << M0->norm << endl;
    delete P;
    delete M0;
}

void NuOsc::cal_P(const FieldVar *inField, Pol *inP)
{
    // Calculate the polarization(\vec {P}) of nu and anu and respective |\vec{P}|
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

void NuOsc::cal_Mn(M *inMn, const Pol *inP, unsigned int n)
{
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

void NuOsc::    survival_prob(const FieldVar *inF, const FieldVar *inF0, std::ofstream &surv_prob_ofstream,
                                   const int time)
{
    double num_Pee = 0;
    double num_Pbee = 0;

    double dnom_Pee = 0;
    double dnom_Pbee = 0;

    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            num_Pee += inF->ee[idx(i, j)] * dz * dv;
            num_Pbee += inF->bee[idx(i, j)] * dz * dv;
            dnom_Pee += inF0->ee[idx(i, j)] * dz * dv;
            dnom_Pbee += inF0->bee[idx(i, j)] * dz * dv;
        }
    }
    surv_prob_ofstream << time << "\t"
                       << std::scientific
                       << num_Pee / dnom_Pee << "\t"
                       << num_Pbee / dnom_Pbee << std::endl;
}

/*---------------------------------------------------------------------------*/

void NuOsc::dom_averaged_survival_prob(const FieldVar *inF, const FieldVar *inF0, const int time)
{
    std::ofstream av_spv_ofstream;
    std::string av_spv_fname = "dom_avrgd_surv_prob_" + std::to_string(time) + "_.dat";
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
        double P_ee = 0;
        double P_bee = 0;

        double P_ee0 = 0;
        double P_bee0 = 0;

        for (int j = 0; j < nz; j++)
        {
            P_ee += inF->ee[idx(i, j)] * dz;
            P_bee += inF->bee[idx(i, j)] * dz;

            P_ee0 += inF0->ee[idx(i, j)] * dz;
            P_bee0 += inF0->bee[idx(i, j)] * dz;
        }
        av_spv_ofstream << std::scientific << vz[i] << "\t"
                        << (P_ee / P_ee0) << "\t"
                        << (P_bee / P_bee0) << std::endl;
    }
    av_spv_ofstream.close();
}

/*---------------------------------------------------------------------------*/

void NuOsc::copy_state(const FieldVar *inF, FieldVar *F)
{
    for (int i = 0; i < nvz; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            F->ee[idx(i, j)] = inF->ee[idx(i, j)];
            F->xx[idx(i, j)] = inF->xx[idx(i, j)];
            F->ex_re[idx(i, j)] = inF->ex_re[idx(i, j)];
            F->ex_im[idx(i, j)] = inF->ex_im[idx(i, j)];

            F->bee[idx(i, j)] = inF->bee[idx(i, j)];
            F->bxx[idx(i, j)] = inF->bxx[idx(i, j)];
            F->bex_re[idx(i, j)] = inF->bex_re[idx(i, j)];
            F->bex_im[idx(i, j)] = inF->bex_im[idx(i, j)];
        }
    }
}
/*---------------------------------------------------------------------------*/
