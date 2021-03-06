#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <execution>
#include <numbers>

#include <fftw3.h>

#include <grid/algorithm.hpp>
#include <grid/core.hpp>

#include <wavefield/tags.hpp>

namespace wavefield
{

using namespace std::literals::complex_literals;

template <class grid_vector, class length_t>
void fresnel_fft_impl(
    grid_vector& real, grid_vector& reciprocal,
    length_t wavelength, fftw_plan& plan)
{
    length_t distance
        = (real.range(0).max - real.range(0).min)
          * (reciprocal.range(1).max - reciprocal.range(1).min)
          / (real.range(0).N * wavelength);

    for (auto [xi, eta] : real.lines()) {
        real.at(xi, eta) *= std::polar(1.0, std::numbers::pi * (xi * xi + eta * eta) / (wavelength * distance));
    }

    fftw_execute(plan);

    Grid::fftshift(reciprocal);

    auto coef = real.range(0).cell_size * real.range(1).cell_size
                / (1.0i * wavelength * distance)
                / (double)(real.range(0).N * real.range(1).N);

    for (auto [x, y] : reciprocal.lines()) {
        reciprocal.at(x, y) *= coef * std::polar(1.0, std::numbers::pi * (distance + (x * x + y * y) / wavelength) / distance);
    }
}

template <class grid_vector, class length_t>
class FresnelFFTDiffraction
{
    grid_vector& real;
    grid_vector& reciprocal;

    length_t wavelength;

    fftw_plan plan;

public:
    FresnelFFTDiffraction(
        grid_vector& real,
        grid_vector& reciprocal,
        length_t wavelength) : real(real), reciprocal(reciprocal), wavelength(wavelength)
    {
        plan = fftw_plan_dft_2d(real.range(0).N, real.range(1).N,
            reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_FORWARD, FFTW_MEASURE);
    }

    void reconfigure()
    {
        plan = fftw_plan_dft_2d(real.range(0).N, real.range(1).N,
            reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_FORWARD, FFTW_MEASURE);
    }

    void propagate()
    {
        fresnel_fft_impl(real, reciprocal, wavelength, plan);
    }

    ~FresnelFFTDiffraction()
    {
        fftw_destroy_plan(plan);
    }
};

template <class grid_vector, class length_t>
void propagate(fresnel_fft_diffraction_tag,
    grid_vector& real, grid_vector& reciprocal,
    length_t wavelength)
{
    fftw_plan plan = fftw_plan_dft_2d(real.range(0).N, real.range(1).N,
        reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_FORWARD, FFTW_ESTIMATE);

    fresnel_fft_impl(real, reciprocal, wavelength, plan);
    fftw_destroy_plan(plan);
}

}  // namespace wavefield
