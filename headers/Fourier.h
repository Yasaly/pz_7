#pragma once
#ifndef Discrete_Fourier_Transformation_h
#define Discrete_Fourier_Transformation_h

#include <vector>
#include <complex>
#include "consts.h"

using namespace std;

namespace Com_Methods
{
    class Fourier
    {
    public:
        void FFT(const vector<complex<double>>& Data, vector<complex<double>>& Result);
        void IFFT(const vector<complex<double>>& Data, vector<complex<double>>& Result);
        void Convolution(const vector<complex<double>>& Vec1,
                         const vector<complex<double>>& Vec2,
                         vector<complex<double>>& Result);
    };
}
#endif

