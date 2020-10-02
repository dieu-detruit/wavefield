#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <execution>
#include <numbers>

#include <fftw3.h>

#include <grid/algorithm.hpp>
#include <grid/bundle.hpp>
#include <grid/core.hpp>

#include <wavefield/tags.hpp>

namespace wavefield
{

using namespace std::literals::complex_literals;

template <class amp_complex_type, class length_t>
class AngularSpectrumDiffraction
{
    using real_grid = Grid::GridVector<amp_complex_type, length_t, 2>;
    using fft_grid = Grid::GridVector<amp_complex_type, decltype(1.0 / length_t{}), 2>;

    real_grid& real;
    real_grid real_expanded;
    fft_grid real_fft;
    real_grid reciprocal_expanded;
    real_grid& reciprocal;

    length_t wavelength;

    fftw_plan plan_fft;
    fftw_plan plan_ifft;

public:
    AngularSpectrumDiffraction(
        real_grid& real,
        real_grid& reciprocal,
        length_t wavelength)
        : real(real),
          real_expanded{real.range(0) * 2, real.range(1) * 2},
          real_fft{{-1.0 / real.range(0).cell_size, 1.0 / real.range(0).cell_size, 2 * real.range(0).N},
              {-1.0 / real.range(1).cell_size, 1.0 / real.range(1).cell_size, 2 * real.range(1).N}},
          reciprocal_expanded{reciprocal.range(0) * 2, reciprocal.range(1) * 2},
          reciprocal(reciprocal),
          wavelength(wavelength)
    {
        plan_fft = fftw_plan_dft_2d(real_expanded.range(0).N, real_expanded.range(1).N,
            reinterpret_cast<fftw_complex*>(real_expanded.data()), reinterpret_cast<fftw_complex*>(real_fft.data()), FFTW_FORWARD, FFTW_MEASURE);

        plan_ifft = fftw_plan_dft_2d(real_fft.range(0).N, real_fft.range(1).N,
            reinterpret_cast<fftw_complex*>(real_fft.data()), reinterpret_cast<fftw_complex*>(reciprocal_expanded.data()), FFTW_BACKWARD, FFTW_MEASURE);
    }

    void propagate(bool forward = true)
    {
        real_expanded.fill(0.0 * *real.begin());
        for (auto [x, y, r] : Grid::zip(real.lines(), real)) {
            real_expanded.at(x, y) = r;
        }

        // 実空間波面をフーリエ変換
        fftw_execute(plan_fft);

        Grid::fftshift(real_fft);

        // 畳み込み(フーリエ変換先)
        {
            length_t distance
                = (forward ? 1.0 : -1.0)
                  * (real.range(0).max - real.range(0).min)
                  * (reciprocal.range(1).max - reciprocal.range(1).min)
                  / (real.range(0).N * wavelength);
            auto coef = 1.0 / (real.range(0).N * real.range(1).N);

            auto urange = real_fft.range(0) / 2;
            auto vrange = real_fft.range(1) / 2;
            for (auto [u, v] : real_fft.lines()) {
                auto w_sq = 1.0 / (wavelength * wavelength) - u * u - v * v;
                if (w_sq > 0.0 * w_sq) {
                    real_fft.at(u, v) *= coef * std::polar(1.0, 2.0 * std::numbers::pi * std::sqrt(w_sq) * distance);
                }
            }
        }

        // 逆フーリエ変換して逆空間波面を得る
        fftw_execute(plan_ifft);

        for (auto [x, y, r] : Grid::zip(reciprocal.lines(), reciprocal)) {
            r = reciprocal_expanded.at(x, y);
        }
    }

    void back_propagate()
    {
        propagate(false);
    }

    ~AngularSpectrumDiffraction()
    {
        fftw_destroy_plan(plan_fft);
        fftw_destroy_plan(plan_ifft);
    }
};  // namespace wavefield

template <class grid_vector, class length_t>
AngularSpectrumDiffraction(grid_vector&, grid_vector&, length_t) -> AngularSpectrumDiffraction<typename grid_vector::value_t, length_t>;

template <class grid_vector, class length_t>
void propagate(angular_spectrum_diffraction_tag,
    grid_vector& real, grid_vector& reciprocal,
    length_t wavelength)
{
}

}  // namespace wavefield
