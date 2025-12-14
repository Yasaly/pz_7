#include "Fourier.h"
#include <stdexcept>
#include <cmath>

using namespace std;

namespace Com_Methods
{
    static inline bool IsPowerOfTwo(size_t n)
    {
        return n && ((n & (n - 1)) == 0);
    }

    static inline size_t ReverseBits(size_t x, unsigned bits)
    {
        size_t y = 0;
        for (unsigned i = 0; i < bits; ++i)
        {
            y = (y << 1) | (x & 1u);
            x >>= 1u;
        }
        return y;
    }
    void Fourier::FFT(const vector<complex<double>>& Data,
                                              vector<complex<double>>& Result)
    {
        const size_t N = Data.size();
        if (!IsPowerOfTwo(N))
            throw runtime_error("FFT error: N must be power of two");

        Result.assign(N, {0.0, 0.0});
        unsigned bits = 0;
        while ((1u << bits) < N) ++bits;

        for (size_t i = 0; i < N; ++i)
        {
            Result[ReverseBits(i, bits)] = Data[i];
        }
        for (size_t len = 2; len <= N; len <<= 1)
        {
            const double ang = -2.0 * PI / static_cast<double>(len);
            const complex<double> wlen(cos(ang), sin(ang));

            for (size_t i = 0; i < N; i += len)
            {
                complex<double> w(1.0, 0.0);
                const size_t half = len >> 1;

                for (size_t j = 0; j < half; ++j)
                {
                    const complex<double> u = Result[i + j];
                    const complex<double> v = Result[i + j + half] * w;
                    Result[i + j] = u + v;
                    Result[i + j + half] = u - v;
                    w *= wlen;
                }
            }
        }
    }

    void Fourier::IFFT(const vector<complex<double>>& Data,
                                               vector<complex<double>>& Result)
    {
        const size_t N = Data.size();
        if (!IsPowerOfTwo(N))
            throw runtime_error("IFFT error: N must be power of two");
        vector<complex<double>> tmp(N);
        for (size_t i = 0; i < N; ++i) tmp[i] = conj(Data[i]);

        FFT(tmp, Result);

        const double invN = 1.0 / static_cast<double>(N);
        for (size_t i = 0; i < N; ++i) Result[i] = conj(Result[i]) * invN;
    }
    void Fourier::Convolution(const vector<complex<double>>& Vec1,
                                                      const vector<complex<double>>& Vec2,
                                                      vector<complex<double>>& Result)
    {
        const size_t N = Vec1.size();
        if (Vec2.size() != N) throw runtime_error("Convolution error: sizes must match");
        if (!IsPowerOfTwo(N)) throw runtime_error("Convolution error: N must be power of two");

        vector<complex<double>> A, B;
        FFT(Vec1, A);
        FFT(Vec2, B);

        for (size_t i = 0; i < N; ++i) A[i] *= B[i];

        IFFT(A, Result);
    }
}
