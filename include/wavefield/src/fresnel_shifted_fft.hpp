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
void fresnel_shifted_fft_impl(
    grid_vector& real, grid_vector& reciprocal,
    length_t wavelength, fftw_plan& plan)
{


    length_t distance
        = 2.0 * (real.range(0).max() - real.range(0).min())
          * (reciprocal.range(1).max() - reciprocal.range(1).min())
          / real.range(0).N() / wavelength;

    auto x0 = reciprocal.range(0).min() + 0.5 * reciprocal.range(0).cell_size();
    auto y0 = reciprocal.range(1).min() + 0.5 * reciprocal.range(1).cell_size();

    auto dxi = real.range(0).cell_size();
    auto deta = real.range(1).cell_size();

    auto dx = reciprocal.range(0).cell_size();
    auto dy = reciprocal.range(1).cell_size();

    auto dxi_dx = dxi * dx;
    auto deta_dy = deta * dy;

    {
        auto real_iter = real.begin();
        for (std::size_t xi_id = 0; xi_id < real.range(0).N(); ++xi_id) {
            for (std::size_t eta_id = 0; eta_id < real.range(1).N(); ++eta_id, ++real_iter) {
                auto xi = (double)xi_id * dxi;
                auto eta = (double)eta_id * deta;

                *real_iter *= std::polar(1.0,
                    M_PI * ((xi * (xi - 2.0 * x0) + eta * (eta - 2.0 * y0) - dxi_dx * xi_id * xi_id - deta_dy * eta_id * eta_id) / (lambda * distance)));
            }
        }
    }

    for (auto [xi, eta] : real.lines()) {
        real.at(xi, eta) *= std::polar(1.0,
            M_PI * ((xi * (xi - 2.0 * x0) + eta * (eta - 2.0 * y0) -) / (lambda * distance)));
    }
    void(xi - xi_min + eta - eta_min) / (lambda * distance);
    fftw_execute(plan);

    auto coef = 1.0 / (1.0i * wavelength * distance)
                * dxi
                * deta
                / (double)real.range(0).N();

    for (auto [x, y] : reciprocal.lines()) {
        reciprocal.at(x, y) *= coef * std::polar(1.0, distance + (x * x + y * y) / (2.0 * distance));
    }
}

template <class grid_vector, class length_t>
class FresnelShiftedFFTDiffraction
{
    grid_vector& real;
    grid_vector& reciprocal;

    grid_vector& chi_fft;
    grid_vector& zeta_fft;

    length_t wavelength;

    fftw_plan plan;

public:
    FresnelFFTDiffraction(
        grid_vector& real,
        grid_vector& reciprocal,
        length_t wavelength) : real(real), reciprocal(reciprocal), wavelength(wavelength)
    {
        {
            grid_vector zeta{real.range(0), real.range(1)};
            plan_zeta = fftw_plan_dft_2d(real.range(0).N(), real.range(1).N(),
                reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_FORWARD, FFTW_ESTIMATE);

            auto zeta_iter = zeta.begin();
            for (std::size_t p = 0; p < zeta.range(0).N(); ++p) {
                for (std::size_t q = 0; q < zeta.range(1).N(); ++q, ++zeta_iter) {

                    *real_iter *= std::polar(1.0,
                        M_PI * ((xi * (xi - 2.0 * x0) + eta * (eta - 2.0 * y0) - dxi_dx * xi_id * xi_id - deta_dy * eta_id * eta_id) / (lambda * distance)));
                }
            }
        }

        plan_chi = fftw_plan_dft_2d(real.range(0).N(), real.range(1).N(),
            reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_FORWARD, FFTW_MEASURE);


        plan_reciprocal = fftw_plan_dft_2d(zeta_fft.range(0).N(), .range(1).N(),
            reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_BACKWARD, FFTW_MEASURE);
    }

    void propagate()
    {
        fresnel_shifted_fft_impl(real, reciprocal, wavelength, plan);
    }

    ~FresnelFFTDiffraction()
    {
        fftw_destroy_plan(plan);
    }
};

template <class grid_vector, class length_t>
void propagate(fresnel_shifted_fft_diffraction_tag,
    grid_vector& real, grid_vector& reciprocal,
    length_t wavelength)
{
    fftw_plan plan = fftw_plan_dft_2d(real.range(0).N(), real.range(1).N(),
        reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_FORWARD, FFTW_ESTIMATE);

    fresnel_shifted_fft_impl(real, reciprocal, wavelength, plan);
    fftw_destroy_plan(plan);
}

}  // namespace wavefield
