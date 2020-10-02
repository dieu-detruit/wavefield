#include <complex>
#include <iostream>

#include <wavefield/core.hpp>

int main()
{
    constexpr double lambda = 4.13e-9;
    constexpr double d = 0.001;
    constexpr double exit_length = lambda * d / 1.0 / 256.0;

    Grid::GridVector<std::complex<double>, double, 2> exit{{-exit_length, exit_length, 256}, {-exit_length, exit_length, 256}};
    Grid::GridVector<std::complex<double>, double, 2> detector{{-1.0, 1.0, 256}, {-1.0, 1.0, 256}};

    wavefield::FraunhoferDiffraction real_to_reciprocal{exit, detector, lambda};

    auto range = Grid::arange(-0.01 * exit_length, 0.01 * exit_length, 2.0 * exit_length / 256.0);
    exit.fill(0.0);
    for (auto [x, y] : Grid::prod(range, range)) {
        exit.at(x, y) = 1.0;
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
