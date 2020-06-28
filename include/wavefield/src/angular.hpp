#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <execution>

#include <fftw3.h>

#include <grid/algorithm.hpp>
#include <grid/core.hpp>

#include <wavefield/tags.hpp>

namespace wavefield
{

using namespace std::literals::complex_literals;

template <class grid_vector, class length_t>
class AngularSpectrumFFTDiffraction
{
    grid_vector& real;
    grid_vector& real_expand;
    grid_vector& real_fft;
    grid_vector& reciprocal;

    length_t wavelength;

    fftw_plan plan_real_fft;

public:
    AngularSpectrumDiffraction(
        grid_vector& real,
        grid_vector& reciprocal,
        length_t wavelength) : real(real),
                               real_expand{{real.range(0).min, real.range(0).max, real.range(0).N}, real.range(1)},
                               real_fft{{real.range(0).min, real.range(0).min}, real.range(1)},
                               reciprocal(reciprocal), wavelength(wavelength)
    {
        plan_real_fft = fftw_plan_dft_2d(real.range(0).N, real.range(1).N,
            reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(real_fft.data()), FFTW_FORWARD, FFTW_MEASURE);
    }

    void propagate_impl()
    {
        // 実空間波面をフーリエ変換
        fftw_execute(plan_real_fft);
        Grid::fftshift(plan_real_fft);

        auto Lx = real.range(0).max - real.range(0).min;
        auto Ly = real.range(1).max - real.range(1).min;

        // 畳み込み(フーリエ変換先)
        {
            auto real_fft_iter = real_fft.begin();
            for (std::size_t i = 0; i < real.range(0).N; i++) {
                for (std::size_t j = 0; j < real.range(1).N; j++, ++real_fft_iter) {
                    auto kx = (i - real.range(0).N / 2.0) / Lx;
                    auto ky = (j - real.range(1).N / 2.0) / Ly;
                    auto kz_normalized = std::sqrt(1.0 - kx * kx - ky * ky);

                    *real_fft_iter *= exp(1.0i * 2.0 * pi * kz * z);
                }
            }
        }

        Grid::fftshift(plan_real_fft);

        length_t distance
            = 2.0 * (real.range(0).max - real.range(0).min)
              * (reciprocal.range(1).max - reciprocal.range(1).min)
              / real.range(0).N / wavelength;

        Grid::fftshift(reciprocal);

        auto coef = real.range(0).cell_size
                    * real.range(1).cell_size
                    / (double)real.range(0).N;

        for (auto [x, y] : reciprocal.lines()) {
            reciprocal.at(x, y) *= coef * std::polar(1.0, k * (distance + (x * x + y * y) / (2.0 * distance)));
        }
    }

    void propagate()
    {
        angular_impl(real, reciprocal, wavelength, plan);
    }

    ~FresnelFFTDiffraction()
    {
        fftw_destroy_plan(plan);
    }
};

template <class grid_vector, class length_t>
void propagate(angular_diffraction_tag,
    grid_vector& real, grid_vector& reciprocal,
    length_t wavelength)
{
    fftw_plan plan = fftw_plan_dft_2d(real.range(0).N, real.range(1).N,
        reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_FORWARD, FFTW_ESTIMATE);

    angular_impl(real, reciprocal, wavelength, plan);
    fftw_destroy_plan(plan);
}

}  // namespace wavefield
