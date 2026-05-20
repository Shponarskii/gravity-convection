#pragma once

#include <vector>

double d0x(
    const std::vector<std::vector<double>>& F,
    int i,
    int j,
    double hx
);

double d0y(
    const std::vector<std::vector<double>>& F,
    int i,
    int j,
    double hy
);

double upwindX(
    const std::vector<std::vector<double>>& F,
    const std::vector<std::vector<double>>& Vx,
    int i,
    int j,
    double hx
);

double upwindY(
    const std::vector<std::vector<double>>& F,
    const std::vector<std::vector<double>>& Vy,
    int i,
    int j,
    double hy
);

double laplacian(
    const std::vector<std::vector<double>>& F,
    int i,
    int j,
    double hx,
    double hy
);