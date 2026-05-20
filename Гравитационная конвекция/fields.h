#pragma once
#ifndef FIELDS_H
#define FIELDS_H

#include <vector>
using namespace std;

struct Fields {
    int M, N;

    vector<vector<double>> C, T, phi, psi, xi;
    vector<vector<double>> Vx, Vy;
    vector<vector<double>> Ex, Ey;
    vector<vector<double>> ix, iy;

    Fields(int M_, int N_);
    void initialize();
};

#endif