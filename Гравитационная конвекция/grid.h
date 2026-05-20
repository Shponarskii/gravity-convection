#pragma once
#ifndef GRID_H
#define GRID_H

struct Grid {
    int M, N;
    double hx, hy;
    double ht;
    double htau1, htau2;
    double L;

    Grid(int M_, int N_, double L_, double ht_, double htau1_, double htau2_);
};

#endif