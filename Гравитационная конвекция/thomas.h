#pragma once
#ifndef THOMAS_H
#define THOMAS_H

#include <vector>
using namespace std;

vector<double> thomasSolver(
    const vector<double>& A,
    const vector<double>& B,
    const vector<double>& D,
    const vector<double>& R
);

#endif