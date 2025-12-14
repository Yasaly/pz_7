#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <iomanip>
#include <cmath>

#include "Wavelet_Analysis.h"

static void SaveSignalCSV(const std::string& path,
                          const std::vector<std::complex<double>>& z,
                          const std::vector<std::complex<double>>& P,
                          const std::vector<std::complex<double>>& Q,
                          const std::vector<std::complex<double>>& rec)
{
    std::ofstream out(path);
    out << "n,z_re,z_im,P_re,P_im,Q_re,Q_im,rec_re,rec_im\n";
    const size_t N = z.size();
    for (size_t i = 0; i < N; ++i)
    {
        out << i << ","
            << z[i].real() << "," << z[i].imag() << ","
            << P[i].real() << "," << P[i].imag() << ","
            << Q[i].real() << "," << Q[i].imag() << ","
            << rec[i].real() << "," << rec[i].imag() << "\n";
    }
}

static void SaveCoefCSV(const std::string& path,
                        int stage,
                        const std::vector<std::complex<double>>& cPsi,
                        const std::vector<std::complex<double>>& cFi)
{
    std::ofstream out(path);
    out << "k,pos,psi_re,psi_im,psi_abs,fi_re,fi_im,fi_abs\n";
    const size_t M = cPsi.size();
    const int step = static_cast<int>(std::pow(2, stage));
    for (size_t k = 0; k < M; ++k)
    {
        const int pos = step * static_cast<int>(k);
        const double psiAbs = std::abs(cPsi[k]);
        const double fiAbs  = std::abs(cFi[k]);
        out << k << "," << pos << ","
            << cPsi[k].real() << "," << cPsi[k].imag() << "," << psiAbs << ","
            << cFi[k].real()  << "," << cFi[k].imag()  << "," << fiAbs  << "\n";
    }
}

static std::string BasisName(Com_Methods::Wavelet_Analysis::Basis_Type t)
{
    using BT = Com_Methods::Wavelet_Analysis::Basis_Type;
    switch (t)
    {
        case BT::Haar: return "Haar";
        case BT::Complex_Shannon: return "Shannon";
        case BT::Daubechies_D6: return "D6";
        default: return "Unknown";
    }
}

//4
static void SavePartialPCSV(const std::string& path,
                            const std::vector<std::complex<double>>& z,
                            const std::vector<std::complex<double>>& P)
{
    std::ofstream out(path);
    out << "n,z_re,z_im,P_re,P_im\n";
    const size_t N = z.size();
    for (size_t i = 0; i < N; ++i)
    {
        out << i << ","
            << z[i].real() << "," << z[i].imag() << ","
            << P[i].real() << "," << P[i].imag() << "\n";
    }
}

static void SaveFilteredCSV(const std::string& path,
                            const std::vector<std::complex<double>>& z,
                            const std::vector<std::complex<double>>& zf)
{
    std::ofstream out(path);
    out << "n,z_re,z_im,zf_re,zf_im\n";
    for (size_t i = 0; i < z.size(); ++i)
    {
        out << i << ","
            << z[i].real() << "," << z[i].imag() << ","
            << zf[i].real() << "," << zf[i].imag() << "\n";
    }
}

//6
static void GenerateSignal6(int N, double A, double B, int w1, int w2, double phi,
                            std::vector<std::complex<double>>& z)
{
    z.assign(N, std::complex<double>(0.0, 0.0));

    for (int j = 0; j < N; ++j)
    {
        double val = 0.0;
        val += A * std::cos(2.0 * PI * w1 * j / double(N) + phi);
        val += B * std::cos(2.0 * PI * w2 * j / double(N));
        z[j] = std::complex<double>(val, 0.0);
    }
}
static void SaveSignalCSV(const std::string& path,
                          const std::vector<std::complex<double>>& z)
{
    std::ofstream out(path);
    out << "n,z_re,z_im\n";
    for (size_t i = 0; i < z.size(); ++i)
        out << i << "," << z[i].real() << "," << z[i].imag() << "\n";
}


int main()
{
    const int N = 512;
    const int Stages = 4;

    const double A = 3.00;
    const double B = 0.26;
    const double w2 = 193.0;

    std::vector<std::complex<double>> z(N, {0.0, 0.0});

    for (int j = 0; j < N; ++j)
    {
        const bool seg2 = (j >= N/4 && j <= N/2);
        const bool seg4 = (j > 3*N/4 && j < N); // (3N/4, N)

        if (seg2 || seg4)
        {
            const double val = A + B * std::cos(2.0 * PI * w2 * j / N);
            z[j] = {val, 0.0};
        }
    }

    //3,4
    using BT = Com_Methods::Wavelet_Analysis::Basis_Type;
    std::vector<BT> bases = { BT::Haar, BT::Complex_Shannon, BT::Daubechies_D6 };

    for (auto basis : bases)
    {
        Com_Methods::Wavelet_Analysis wa(N, basis);

        for (int stage = 1; stage <= Stages; ++stage)
        {
            std::vector<std::complex<double>> cPsi, cFi;
            wa.Analysis_Phase(stage, z, cPsi, cFi);

            std::vector<std::complex<double>> P, Q, rec;
            wa.Synthesis_Phase(stage, cPsi, cFi, P, Q, rec);

            const std::string bname = BasisName(basis);
            if (stage <= 3)
            {
                SavePartialPCSV(bname + "_P_minus_" + std::to_string(stage) + ".csv", z, P);
            }

            SaveCoefCSV(bname + "_stage" + std::to_string(stage) + "_coef.csv", stage, cPsi, cFi);
            SaveSignalCSV(bname + "_stage" + std::to_string(stage) + "_recovery.csv", z, P, Q, rec);

            //5
            if (stage == 2)
            {
                auto cPsi_filt = cPsi;
                auto cFi_filt  = cFi;

                for (auto& v : cPsi_filt) v = std::complex<double>(0.0, 0.0);

                std::vector<std::complex<double>> P_f, Q_f, rec_f;
                wa.Synthesis_Phase(stage, cPsi_filt, cFi_filt, P_f, Q_f, rec_f);

                SaveFilteredCSV(bname + "_stage2_filter_P_minus_1.csv", z, rec_f);
            }

            std::cout << bname << " stage " << stage
                      << ": coef=" << cPsi.size()
                      << " saved.\n";
        }
    }

    std::cout << "Done.\n";
    return 0;
}