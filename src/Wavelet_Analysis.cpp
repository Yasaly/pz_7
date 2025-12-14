#include "Wavelet_Analysis.h"

#include <cmath>
#include <stdexcept>
#include <complex>
#include <vector>

using namespace std;

namespace Com_Methods
{
	Wavelet_Analysis::Wavelet_Analysis(int Count_Data, Basis_Type Type)
	{
		const int N = Count_Data;
		if (N <= 0) throw runtime_error("Wavelet_Analysis: N must be positive");

		U.assign(N, complex<double>(0.0, 0.0));
		V.assign(N, complex<double>(0.0, 0.0));

		switch (Type)
		{
			case Basis_Type::Complex_Shannon:
			{
				if (N % 4 != 0)
					throw runtime_error("Complex Shannon Basis: N % 4 != 0");

				U[0] = V[0] = 1.0 / sqrt(2.0);

				for (int i = 1; i < N; ++i)
				{
					const double denom = sin(PI * i / N);
					if (abs(denom) < EPS)
						throw runtime_error("Complex Shannon Basis: division by zero in sin(PI*i/N)");

					const double a = sqrt(2.0) / N * cos(PI * i / N) * sin(PI * i / 2.0) / denom;
					const double b = -sqrt(2.0) / N * sin(PI * i / N) * sin(PI * i / 2.0) / denom;

					const complex<double> val(a, b);
					U[i] = val;
					V[i] = ((i % 2) ? -1.0 : 1.0) * val;
				}
				break;
			}

			case Basis_Type::Haar:
			{
				if (N % 2 != 0)
					throw runtime_error("Haar Basis: N must be even");

				const double s = 1.0 / sqrt(2.0);
				U.assign(N, complex<double>(0.0, 0.0));
				V.assign(N, complex<double>(0.0, 0.0));

				U[0] = s; U[1] = s;
				V[0] = s; V[1] = -s;
				break;
			}

			case Basis_Type::Daubechies_D6:
			{
				if (N <= 6)
					throw runtime_error("D6 Basis: N must be > 6");
				if (N % 2 != 0)
					throw runtime_error("D6 Basis: N must be even");

				const double rt10 = sqrt(10.0);
				const double a = 1.0 - rt10;
				const double b = 1.0 + rt10;
				const double c = sqrt(5.0 + 2.0 * rt10);
				const double k = sqrt(2.0) / 32.0;

				U.assign(N, complex<double>(0.0, 0.0));
				V.assign(N, complex<double>(0.0, 0.0));

				U[0] = k * (b + c);
				U[1] = k * (2.0 * a + 3.0 * b + 3.0 * c);
				U[2] = k * (6.0 * a + 4.0 * b + 2.0 * c);
				U[3] = k * (6.0 * a + 4.0 * b - 2.0 * c);
				U[4] = k * (2.0 * a + 3.0 * b - 3.0 * c);
				U[5] = k * (b - c);

				V[0] = -U[1];
				V[1] =  U[0];

				V[N - 4] = -U[5];
				V[N - 3] =  U[4];
				V[N - 2] = -U[3];
				V[N - 1] =  U[2];

				break;
			}

			default:
				throw runtime_error("Wavelet_Analysis: unknown Basis_Type");
		}
	}

	void Wavelet_Analysis::Wavelet_Filters_System(int Stages)
	{
		if (Stages <= 0) throw runtime_error("Wavelet_Filters_System: Stages must be positive");

		Operators op;
		const int N = static_cast<int>(U.size());

		vector<vector<complex<double>>> U_Filters(Stages), V_Filters(Stages);
		U_Filters[0] = U;
		V_Filters[0] = V;

		for (int i = 1; i < Stages; ++i)
		{
			const int div = static_cast<int>(pow(2.0, i));
			if (N % div != 0)
				throw runtime_error("Wavelet_Filters_System: N must be divisible by 2^i for all stages");

			const int Count_Elements = N / div;
			U_Filters[i].assign(Count_Elements, complex<double>(0.0, 0.0));
			V_Filters[i].assign(Count_Elements, complex<double>(0.0, 0.0));

			for (int n = 0; n < Count_Elements; ++n)
			{
				for (int k = 0; k < div; ++k)
				{
					const int idx = n + k * (N / div);
					U_Filters[i][n] += U_Filters[0][idx];
					V_Filters[i][n] += V_Filters[0][idx];
				}
			}
		}

		Fourier dft;

		vector<complex<double>> U_u, U_v;

		f.resize(Stages);
		g.resize(Stages);

		f[0] = V_Filters[0];
		g[0] = U_Filters[0];

		for (int i = 1; i < Stages; ++i)
		{
			op.Upsampling_Operator(i, U_Filters[i], U_u);
			op.Upsampling_Operator(i, V_Filters[i], U_v);

			dft.Convolution(g[i - 1], U_v, f[i]);
			dft.Convolution(g[i - 1], U_u, g[i]);
		}
	}

	void Wavelet_Analysis::Wavelet_Basis(int Stage,
		vector<vector<complex<double>>>& psi,
		vector<vector<complex<double>>>& phi)
	{
		if (Stage <= 0) throw runtime_error("Wavelet_Basis: Stage must be >= 1");

		Operators op;
		const int N = static_cast<int>(U.size());
		const int step = static_cast<int>(pow(2.0, Stage));

		if (N % step != 0)
			throw runtime_error("Wavelet_Basis: N must be divisible by 2^Stage");

		const int Count_Bas_Elements = N / step;

		if (static_cast<int>(g.size()) < Stage)
			Wavelet_Filters_System(Stage);

		psi.resize(Count_Bas_Elements);
		phi.resize(Count_Bas_Elements);

		for (int i = 0; i < Count_Bas_Elements; ++i)
		{
			const int shift = step * i;
			op.Cyclic_Shift(shift, f[Stage - 1], psi[i]);
			op.Cyclic_Shift(shift, g[Stage - 1], phi[i]);
		}
	}

	void Wavelet_Analysis::Analysis_Phase(int Stage,
		const vector<complex<double>>& Data,
		vector<complex<double>>& coef_psi,
		vector<complex<double>>& coef_phi)
	{
		if (static_cast<int>(Data.size()) != static_cast<int>(U.size()))
			throw runtime_error("Analysis_Phase: Data size mismatch");

		Operators op;

		vector<vector<complex<double>>> psi, phi;
		Wavelet_Basis(Stage, psi, phi);

		const int M = static_cast<int>(psi.size());

		coef_psi.clear();
		coef_phi.clear();
		coef_psi.reserve(M);
		coef_phi.reserve(M);

		for (int k = 0; k < M; ++k)
		{
			coef_psi.push_back(op.Dot_Product(Data, psi[k]));
			coef_phi.push_back(op.Dot_Product(Data, phi[k]));
		}
	}

	void Wavelet_Analysis::Synthesis_Phase(int Stage,
		const vector<complex<double>>& coef_psi,
		const vector<complex<double>>& coef_phi,
		vector<complex<double>>& P,
		vector<complex<double>>& Q,
		vector<complex<double>>& Recovery)
	{
		vector<vector<complex<double>>> psi, phi;
		Wavelet_Basis(Stage, psi, phi);

		const int M = static_cast<int>(psi.size());
		const int N = static_cast<int>(U.size());

		if (static_cast<int>(coef_psi.size()) != M || static_cast<int>(coef_phi.size()) != M)
			throw runtime_error("Synthesis_Phase: coefficient size mismatch");

		P.clear(); Q.clear(); Recovery.clear();
		P.reserve(N); Q.reserve(N); Recovery.reserve(N);

		for (int n = 0; n < N; ++n)
		{
			complex<double> Pn(0.0, 0.0);
			complex<double> Qn(0.0, 0.0);

			for (int k = 0; k < M; ++k)
			{
				Pn += coef_phi[k]  * phi[k][n];
				Qn += coef_psi[k] * psi[k][n];
			}

			P.push_back(Pn);
			Q.push_back(Qn);
			Recovery.push_back(Pn + Qn);
		}
	}
}