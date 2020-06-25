#include <complex>

#include <wavefield/core.hpp>

int main()
{

    constexpr double lambda = 4.13e-9;

    Grid::GridVector<std::complex<double>, double, 2> exit{{-1.0, 1.0, 100}, {-1.0, 1.0, 100}};
    Grid::GridVector<std::complex<double>, double, 2> detector{{-1.0, 1.0, 100}, {-1.0, 1.0, 100}};

    wavefield::FresnelFFTDiffraction real_to_reciprocal{exit, detector, lambda};

    real_to_reciprocal.propagate();

    return 0;
}
