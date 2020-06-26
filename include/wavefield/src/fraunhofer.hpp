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
void fraunhofer_impl(
    grid_vector& real, grid_vector& reciprocal,
    length_t wavelength)
{
    Grid::fftshift(reciprocal);
    length_t distance
        = 2.0 * (real.range(0).max() - real.range(0).min())
          * (reciprocal.range(1).max() - reciprocal.range(1).min())
          / real.range(0).N() / wavelength;

    //auto coef = std::exp(1.0i * wavelength * distance)
    /// (1.0i * wavelength * distance)
    //* real.range(0).cell_size()
    //* real.range(1).cell_size()
    /// (double)real.range(0).N();

    //std::for_each(std::execution::par_unseq,
    //reciprocal.begin(), reciprocal.end(), [&coef](auto& f) {
    //f *= coef;
    //});
}

template <class grid_vector, class length_t>
class FraunhoferDiffraction
{
    grid_vector& real;
    grid_vector& reciprocal;

    length_t wavelength;

    fftw_plan plan;

public:
    FraunhoferDiffraction(
        grid_vector& real,
        grid_vector& reciprocal,
        length_t wavelength) : real(real), reciprocal(reciprocal), wavelength(wavelength)
    {
        plan = fftw_plan_dft_2d(real.range(0).N(), real.range(1).N(),
            reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_BACKWARD, FFTW_MEASURE);
    }

    void propagate()
    {
        fftw_execute(plan);
        fraunhofer_impl(real, reciprocal, wavelength);
    }

    ~FraunhoferDiffraction()
    {
        fftw_destroy_plan(plan);
    }
};

template <class grid_vector, class length_t>
void propagate(fraunhofer_diffraction_tag,
    grid_vector& real, grid_vector& reciprocal,
    length_t wavelength)
{
    fftw_plan plan = fftw_plan_dft_2d(real.range(0).N(), real.range(1).N(),
        reinterpret_cast<fftw_complex*>(real.data()), reinterpret_cast<fftw_complex*>(reciprocal.data()), FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fraunhofer_impl(real, reciprocal, wavelength);
    fftw_destroy_plan(plan);
}

}  // namespace wavefield
