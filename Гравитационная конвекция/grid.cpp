#include "grid.h"

Grid::Grid(int M_, int N_, double L_, double ht_, double htau1_, double htau2_) {
    M = M_;
    N = N_;
    L = L_;
    ht = ht_;
    htau1 = htau1_;
    htau2 = htau2_;

    hx = 1.0 / M;
    hy = L / N;
}