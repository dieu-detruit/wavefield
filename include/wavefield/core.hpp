#pragma once

#include <wavefield/tags.hpp>

#include <wavefield/src/angular.hpp>
#include <wavefield/src/fraunhofer.hpp>
#include <wavefield/src/fresnel_fft.hpp>
#include <wavefield/src/fresnel_shifted_fft.hpp>

namespace wavefield
{
inline constexpr bool forward = true;
inline constexpr bool backward = true;
}  // namespace wavefield
