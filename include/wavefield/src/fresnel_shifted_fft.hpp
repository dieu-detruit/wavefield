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

//template <class grid_vector, class length_t>
//class FresnelShiftedFFTDiffraction
//{

//public:
//FresnelShiftedFFTDiffraction(
//grid_vector& real,
//grid_vector& reciprocal,
//length_t wavelength,
//length_t distance) : real(real), reciprocal(reciprocal), wavelength(wavelength), distance(distance)
//{
//}

//void propagate()
//{
//// u(p, q) -> h(p, q)


//// h(p, q) -> chi(i, j)


//// zeta(i, j) * chi(i, j)


//// U(m, n) = iFFT[zeta(i, j) * chi(i, j)]
//}

//~FresnelShiftedFFTDiffraction()
//{
//fftw_destroy_plan(plan);
//}
//};

//template <class grid_vector, class length_t>
//void propagate(fresnel_shifted_fft_diffraction_tag,
//grid_vector& real, grid_vector& reciprocal,
//length_t wavelength)
//{
//}

}  // namespace wavefield
