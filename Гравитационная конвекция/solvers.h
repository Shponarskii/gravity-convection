#pragma once
#ifndef SOLVERS_H
#define SOLVERS_H

#include "fields.h"
#include "params.h"
#include "grid.h"

void solveC_ADI(Fields& F, const Grid& G, const Params& P);
void solvePhi_ADI(Fields& F, const Grid& G, const Params& P);
void computeE(Fields& F, const Grid& G, const Params& P);
void computeIy(Fields& F, const Grid& G, const Params& P);
void restoreIx(Fields& F, const Grid& G, const Params& P);
void solveT_ADI(Fields& F, const Grid& G, const Params& P);
void solveXi_ADI(Fields& F, const Grid& G, const Params& P);
void solvePsi(Fields& F, const Grid& G, const Params& P);
void updateVelocity(Fields& F, const Grid& G, const Params& P);
#endif