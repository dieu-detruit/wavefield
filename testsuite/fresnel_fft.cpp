#include <complex>
#include <iostream>

#include <wavefield/core.hpp>

int main()
{

    constexpr double lambda = 4.13e-9;
    constexpr double d = 0.0001;

    Grid::GridVector<std::complex<double>, double, 2> exit{{-1.0, 1.0, 256}, {-1.0, 1.0, 256}};
    Grid::GridVector<std::complex<double>, double, 2> detector{{-10.0, 10.0, 256}, {-10.0, 10.0, 256}};

    wavefield::FresnelFFTDiffraction real_to_reciprocal{exit, detector, lambda};

    exit.fill(0.0);
    for (auto [x, y] : exit.lines()) {
        auto r = std::hypot(x, y);
        if (0.1 < r and r < 0.15) {
            exit.at(x, y) = 1.0;
        }
    }

    real_to_reciprocal.propagate();
    for (auto& x : detector.line(0)) {
        for (auto& y : detector.line(1)) {
            std::cout << x << ' ' << y << ' '
                      << std::arg(detector.at(x, y)) << ' '
                      << std::norm(detector.at(x, y)) << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}
