#include "Operators.h"
#include <cmath>

using namespace std;

namespace Com_Methods
{
	void Operators::Cyclic_Shift(int k, const vector<complex<double>>& Data, vector<complex<double>>& Result)
	{
		int N = (int)Data.size();
		Result.clear(); Result.resize(N);

		for (int i = 0; i < N; i++)
		{
			if (i - k < 0)
				Result[i] = Data[i - k + N];
			else
				Result[i] = Data[i - k];
		}
	}

	void Operators::Downsampling_Operator(int l, const vector<complex<double>>& Data, vector<complex<double>>& Result)
	{
		int Val = int(pow(2.0, l));
		int N = (int)Data.size() / Val;
		Result.clear();
		Result.resize(N);
		for (int i = 0; i < N; i++)
		{
			Result[i] = Data[i * Val];
		}
	}
	void Operators::Upsampling_Operator(int l, const vector<complex<double>>& Data, vector<complex<double>>& Result)
	{
		int Val = int(pow(2.0, l));
		int N = (int)Data.size() * Val;
		Result.clear();
		Result.resize(N);
		for (int i = 0; i < N; i++)
		{
			if (i % Val == 0) Result[i] = Data[i / Val];
			else Result[i] = 0;
		}
	}
	complex<double> Operators::Dot_Product(const vector<complex<double>>& Vec1, const vector<complex<double>>& Vec2)
	{
		int N = (int)Vec1.size();
		complex<double> Result(0.0, 0.0);

		for (int i = 0; i < N; i++)
			Result += Vec1[i] * conj(Vec2[i]);

		return Result;
	}
}
