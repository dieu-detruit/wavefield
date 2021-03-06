#pragma once

namespace wavefield
{

struct diffraction_tag {
};

// Fraunhofer
struct fraunhofer_diffraction_tag : diffraction_tag {
};
inline constexpr fraunhofer_diffraction_tag fraunhofer;

// Fresnel
struct fresnel_fft_diffraction_tag : diffraction_tag {
};
inline constexpr fresnel_fft_diffraction_tag fresnel_fft;

struct fresnel_convolution_diffraction_tag : diffraction_tag {
};
inline constexpr fresnel_convolution_diffraction_tag fresnel_convolution;

struct fresnel_shifted_fft_diffraction_tag : diffraction_tag {
};
inline constexpr fresnel_shifted_fft_diffraction_tag fresnel_shifted_fft;

// Angular
struct angular_spectrum_diffraction_tag : diffraction_tag {
};
inline constexpr angular_spectrum_diffraction_tag angular;

}  // namespace wavefield
