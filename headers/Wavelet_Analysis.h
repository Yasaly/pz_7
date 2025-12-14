#pragma once

#ifndef Wavelet_Analysis_H
#define Wavelet_Analysis_H

#include <vector>
#include <complex>
#include "Operators.h"
#include "Fourier.h"

using namespace std;

namespace Com_Methods
{
	class Wavelet_Analysis
	{
	private:
		vector<complex<double>> U, V;
		vector<vector<complex<double>>> f, g;

	public:
		enum class Basis_Type
		{
			Haar = 0,
			Complex_Shannon = 1,
			Daubechies_D6 = 2,
		};
		Wavelet_Analysis(int Count_Data, Basis_Type Type);

	private:
		void Wavelet_Filters_System(int Stages);
		void Wavelet_Basis(int Stage,
						   vector<vector<complex<double>>>& psi,
						   vector<vector<complex<double>>>& phi);

	public:
		void Analysis_Phase(int Stage, 
							const vector<complex<double>>& Data,
							vector<complex<double>>& coef_psi,
							vector<complex<double>>& coef_phi);
		void Synthesis_Phase(int Stage,
								const vector<complex<double>>& coef_psi,
								const vector<complex<double>>& coef_phi,
								vector<complex<double>>& P,
								vector<complex<double>>& Q,
								vector<complex<double>>& Recovery);
	};
}

#endif
