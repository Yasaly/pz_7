#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "consts.h"
#include "Wavelet_Analysis.h"

using namespace std;

using cd = complex<double>;

static string BasisName(Com_Methods::Wavelet_Analysis::Basis_Type b)
{
    using BT = Com_Methods::Wavelet_Analysis::Basis_Type;
    switch (b)
    {
        case BT::Haar:           return "Haar";
        case BT::Complex_Shannon:return "Shannon";
        case BT::Daubechies_D6:  return "D6";
        default:                 return "Unknown";
    }
}

static void GenerateSignalP6(int N, double A, double B, int w1, int w2, double phi, vector<cd>& clean, vector<cd>& noise, vector<cd>& noisy)
{
    clean.assign(N, cd(0.0, 0.0));
    noise.assign(N, cd(0.0, 0.0));
    noisy.assign(N, cd(0.0, 0.0));

    for (int j = 0; j < N; ++j)
    {
        const double s = A * cos(2.0 * PI * w1 * j / double(N) + phi);
        const double n = B * cos(2.0 * PI * w2 * j / double(N));

        clean[j] = cd(s, 0.0);
        noise[j] = cd(n, 0.0);
        noisy[j] = cd(s + n, 0.0);
    }
}

static void SaveP6SignalCSV(const string& path, const vector<cd>& clean, const vector<cd>& noise, const vector<cd>& noisy)
{
    ofstream out(path);
    out << "n,clean,noise,noisy\n";
    for (size_t i = 0; i < noisy.size(); ++i)
        out << i << "," << clean[i].real() << "," << noise[i].real() << "," << noisy[i].real() << "\n";
}

static void SaveCoefCSV(const string& path, const vector<cd>& cPsi, const vector<cd>& cPhi)
{
    ofstream out(path);
    out << "k,abs_psi,abs_phi\n";
    const size_t M = min(cPsi.size(), cPhi.size());
    for (size_t k = 0; k < M; ++k)
        out << k << "," << abs(cPsi[k]) << "," << abs(cPhi[k]) << "\n";
}

static void SaveStagePQRecoveryCSV(const string& path, const vector<complex<double>>& z, const vector<complex<double>>& P,
                                   const vector<complex<double>>& Q, const vector<complex<double>>& recovery)
{
    ofstream out(path);
    out << "n,z_re,P_re,Q_re,recovery_re\n";

    for (size_t i = 0; i < z.size(); ++i)
    {
        out << i << ","
            << z[i].real() << ","
            << P[i].real() << ","
            << Q[i].real() << ","
            << recovery[i].real() << "\n";
    }
}


int main()
{
    const int n = 9;
    const int N = 512;

    const double A = 3.0;
    const double B = 0.26;
    const int w1 = 1;
    const int w2 = 193;
    const double phi = 0.0;

    const int Stages = 4;

    vector<cd> clean, noise, noisy;
    GenerateSignalP6(N, A, B, w1, w2, phi, clean, noise, noisy);
    SaveP6SignalCSV("p6_signal.csv", clean, noise, noisy);

    cout << "p6_signal.csv saved. N=" << N
              << " A=" << A << " B=" << B
              << " w1=" << w1 << " w2=" << w2
              << " phi=" << phi << "\n";

    using WA = Com_Methods::Wavelet_Analysis;
    vector<WA::Basis_Type> bases = {
        WA::Basis_Type::Haar,
        WA::Basis_Type::Complex_Shannon,
        WA::Basis_Type::Daubechies_D6
    };

    for (auto basis : bases)
    {
        const string bname = BasisName(basis);
        WA wa(N, basis);

        for (int stage = 1; stage <= Stages; ++stage)
        {
            vector<cd> cPsi, cPhi;
            wa.Analysis_Phase(stage, noisy, cPsi, cPhi);

            SaveCoefCSV("p6_" + bname + "_stage" + to_string(stage) + "_coef.csv", cPsi, cPhi);

            vector<cd> P, Q, rec;
            wa.Synthesis_Phase(stage, cPsi, cPhi, P, Q, rec);

            SaveStagePQRecoveryCSV(
                "p6_" + bname + "_stage" + to_string(stage) + "_PQ_recovery.csv",
                noisy, P, Q, rec
            );

            cout << bname << " stage=" << stage
                      << " saved: coef + PQ_filter\n";
        }
    }

    cout << "Done.\n";
    return 0;
}
