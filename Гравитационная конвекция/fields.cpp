#include "fields.h"
#include "params.h"

Params P;
void Fields::initialize() {

    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= N; j++) {

            double x = double(i) / M;

            C[i][j] = 1.0;
            T[i][j] = 1.0;

            Vx[i][j] = 0.0;
            if (i == 0 or i == M)
                Vy[i][j] = 0.0;
            else Vy[i][j] = 1.0;

            phi[i][j] = P.d_phi * x;
            psi[i][j] = 1.0 - x;

            xi[i][j] = 0.0;

            Ex[i][j] = -P.d_phi;
            Ey[i][j] = 0.0;
        }
    }
}

Fields::Fields(int M_, int N_) {
    M = M_;
    N = N_;

    C.assign(M + 1, vector<double>(N + 1, 0.0));
    T.assign(M + 1, vector<double>(N + 1, 0.0));
    phi.assign(M + 1, vector<double>(N + 1, 0.0));
    psi.assign(M + 1, vector<double>(N + 1, 0.0));
    xi.assign(M + 1, vector<double>(N + 1, 0.0));

    Vx.assign(M + 1, vector<double>(N + 1, 0.0));
    Vy.assign(M + 1, vector<double>(N + 1, 0.0));

    Ex.assign(M + 1, vector<double>(N + 1, 0.0));
    Ey.assign(M + 1, vector<double>(N + 1, 0.0));

    ix.assign(M + 1, vector<double>(N + 1, 0.0));
    iy.assign(M + 1, vector<double>(N + 1, 0.0));
}