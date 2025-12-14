#pragma once
#ifndef Complex_Operators_h
#define Complex_Operators_h

#include <complex>
#include <vector>
#include "consts.h"

using namespace std;

namespace Com_Methods
{
    class Operators
    {
    public:
        void Cyclic_Shift(int k, const vector<complex<double>>& Data, vector<complex<double>>& Result);
        static void Downsampling_Operator(int l, const vector<complex<double>>& Data, vector<complex<double>>& Result);
        void Upsampling_Operator(int l, const vector<complex<double>>& Data, vector<complex<double>>& Result);
        complex<double> Dot_Product(const vector<complex<double>>& Vec1, const vector<complex<double>>& Vec2);
    };
}

#endif

