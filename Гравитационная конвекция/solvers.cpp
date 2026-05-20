#include <iostream>
#include "solvers.h"
#include "thomas.h"
#include "operators.h"
#include <vector>
using namespace std;

void solveC_ADI(Fields& F, const Grid& G, const Params& P)
{
    vector<vector<double>> Cstar = F.C;
    vector<vector<double>> Cnew = F.C;

    for (int j = 1; j < G.N; j++)
    {
        vector<double> A(G.M + 1, 0.0);
        vector<double> B(G.M + 1, 0.0);
        vector<double> D(G.M + 1, 0.0);
        vector<double> R(G.M + 1, 0.0);
        for (int i = 1; i < G.M; i++)
        {
            double Vx = F.Vx[i][j];
            double Vy = F.Vy[i][j];

            if (Vx >= 0.0)
            {
                A[i] = 1.0 / (P.Pe * G.hx * G.hx) + Vx / G.hx;
                B[i] = -2.0 / G.ht - 2.0 / (P.Pe * G.hx * G.hx) - Vx / G.hx;
                D[i] = 1.0 / (P.Pe * G.hx * G.hx);
            }
            else
            {
                A[i] = 1.0 / (P.Pe * G.hx * G.hx);
                B[i] = -2.0 / G.ht - 2.0 / (P.Pe * G.hx * G.hx) + Vx / G.hx;
                D[i] = 1.0 / (P.Pe * G.hx * G.hx) - Vx / G.hx;
            }
            double dCdy = upwindY(F.C, F.Vy, i, j, G.hy);
            double diffusionY = (F.C[i][j + 1] - 2.0 * F.C[i][j] + F.C[i][j - 1]) / (G.hy * G.hy);

            R[i] = -2.0 * F.C[i][j] / G.ht + Vy * dCdy - diffusionY / P.Pe;
        }
        double ix_left = (P.KI1 * F.C[0][j] * F.Ex[0][j]) / (1.0 - P.KI2 * P.tA);
        double ix_right = (P.KI1 * F.C[G.M][j] * F.Ex[G.M][j]) / (1.0 + P.KI2 * P.tK);

        A[0] = 0.0;
        B[0] = 1.0;
        D[0] = -1.0;
        R[0] = -G.hx * P.tA * ix_left;

        A[G.M] = -1.0;
        B[G.M] = 1.0;
        D[G.M] = 0.0;
        R[G.M] = -G.hx * P.tK * ix_right;
        vector<double> sol = thomasSolver(A, B, D, R);

        for (int i = 0; i <= G.M; i++)
            Cstar[i][j] = sol[i];
    }

    for (int i = 1; i < G.M; i++)
    {
        vector<double> A(G.N + 1, 0.0);
        vector<double> B(G.N + 1, 0.0);
        vector<double> D(G.N + 1, 0.0);
        vector<double> R(G.N + 1, 0.0);
        for (int j = 1; j < G.N; j++)
        {
            double Vx = F.Vx[i][j];
            double Vy = F.Vy[i][j];

            if (Vy >= 0.0)
            {
                A[j] = 1.0 / (P.Pe * G.hy * G.hy) + Vy / G.hy;
                B[j] = -2.0 / G.ht - 2.0 / (P.Pe * G.hy * G.hy) - Vy / G.hy;
                D[j] = 1.0 / (P.Pe * G.hy * G.hy);
            }
            else
            {
                A[j] = 1.0 / (P.Pe * G.hy * G.hy);
                B[j] = -2.0 / G.ht - 2.0 / (P.Pe * G.hy * G.hy) + Vy / G.hy;
                D[j] = 1.0 / (P.Pe * G.hy * G.hy) - Vy / G.hy;
            }
            double dCdx = upwindX(Cstar, F.Vx, i, j, G.hx);
            double diffusionX = (Cstar[i + 1][j] - 2.0 * Cstar[i][j] + Cstar[i - 1][j]) / (G.hx * G.hx);

            R[j] = -2.0 * Cstar[i][j] / G.ht + Vx * dCdx - diffusionX / P.Pe;
        }

        A[0] = 0.0;
        B[0] = 1.0;
        D[0] = 0.0;
        R[0] = 1.0;

        A[G.N] = -1.0;
        B[G.N] = 1.0;
        D[G.N] = 0.0;
        R[G.N] = 0.0;
        vector<double> sol = thomasSolver(A, B, D, R);

        for (int j = 0; j <= G.N; j++)
            Cnew[i][j] = sol[j];
    }
    F.C = Cnew;
    for (int j = 0; j <= G.N; ++j)
    {
        double ixL = (P.KI1 * F.C[0][j] * F.Ex[0][j]) / (1.0 - P.KI2 * P.tA);
        F.C[0][j] = F.C[1][j] - G.hx * P.tA * ixL;

        double ixR = (P.KI1 * F.C[G.M][j] * F.Ex[G.M][j]) / (1.0 + P.KI2 * P.tK);
        F.C[G.M][j] = F.C[G.M - 1][j] - G.hx * P.tK * ixR;
    }
}

void solvePhi_ADI(Fields& F, const Grid& G, const Params& P)
{
    vector<vector<double>> phi_star = F.phi;
    vector<vector<double>> phi_new = F.phi;

    double err = 1.0;
    int iter = 0;
    int iterMax = 100;
    double K_phix = 4 * G.hx * G.hx / G.htau2;
    double beta = G.hx / G.hy;
    double K_phiy = 4 * G.hy * G.hy / G.htau2;
    while (err > P.eps && iter < iterMax)
    {
        vector<vector<double>> phi_old = F.phi;
        for (int j = 1; j < G.N; ++j)
        {
            vector<double> A(G.M + 1), B(G.M + 1), D(G.M + 1), R(G.M + 1);
            for (int i = 1; i < G.M; ++i)
            {
                A[i] = F.C[i][j] + F.C[i - 1][j];
                B[i] = -F.C[i+1][j] - 2*F.C[i][j] - F.C[i-1][j] - K_phix;
                D[i] = F.C[i + 1][j] + F.C[i][j];
                R[i] = pow(beta, 2) * ((F.C[i][j] + F.C[i][j - 1]) * (F.phi[i][j] - F.phi[i][j - 1]) -
                    (F.C[i][j + 1] + F.C[i][j]) * (F.phi[i][j + 1] - F.phi[i][j])) - K_phix * F.phi[i][j]
                    - 2 * pow(G.hx, 2) * P.Kphi * laplacian(F.C, i, j, G.hx, G.hy);
            }

            A[0] = 0; B[0] = 1; D[0] = 0; R[0] = 0;
            A[G.M] = 0; B[G.M] = 1; D[G.M] = 0; R[G.M] = P.d_phi;

            auto sol = thomasSolver(A, B, D, R);
            for (int i = 0; i <= G.M; ++i)
                phi_star[i][j] = sol[i];
        }

        for (int i = 1; i < G.M; ++i)
        {
            vector<double> A(G.N + 1), B(G.N + 1), D(G.N + 1), R(G.N + 1);

            for (int j = 1; j < G.N; ++j)
            {
                A[j] = F.C[i][j] + F.C[i][j-1];
                B[j] = -F.C[i][j + 1] - 2 * F.C[i][j] - F.C[i][j - 1] - K_phiy;
                D[j] = F.C[i][j + 1] + F.C[i][j];
                R[j] = ((F.C[i][j] + F.C[i - 1][j]) * (phi_star[i][j] - phi_star[i - 1][j]) -
                    (F.C[i + 1][j] + F.C[i][j]) * (phi_star[i + 1][j] - phi_star[i][j])) / pow(beta, 2)
                    - K_phiy * phi_star[i][j] - 2 * pow(G.hy, 2) * P.Kphi * laplacian(F.C, i, j, G.hx, G.hy);
            }

            A[0] = 0; B[0] = 1; D[0] = 0; R[0] = P.d_phi * i * G.hx;
            A[G.N] = -1; B[G.N] = 1; D[G.N] = 0; R[G.N] = 0;

            auto sol = thomasSolver(A, B, D, R);

            for (int j = 0; j <= G.N; ++j)
                phi_new[i][j] = sol[j];
        }

        for (int j = 0; j <= G.N; ++j)
        {
            phi_new[0][j] = 0.0;
            phi_new[G.M][j] = P.d_phi;
        }

        err = 0.0;

        for (int i = 0; i <= G.M; ++i)
            for (int j = 0; j <= G.N; ++j)
                err = max(err, fabs(phi_new[i][j] - phi_old[i][j]));

        cout << "Iter = " << iter << " errPhi = " << err << endl;
        F.phi = phi_new;
        iter++;
    }
}

void computeE(Fields& F, const Grid& G, const Params& P)
{
    for (int i = 1; i < G.M; ++i)
    {
        for (int j = 1; j < G.N; ++j)
        {
            F.Ex[i][j] = -d0x(F.phi, i, j, G.hx);
            F.Ey[i][j] = -d0y(F.phi, i, j, G.hy);
        }
    }

    for (int j = 0; j <= G.N; ++j)
    {
        F.Ex[0][j] = (F.phi[0][j] - F.phi[1][j]) / G.hx;
        F.Ey[0][j] = 0.0;

        F.Ex[G.M][j] = (F.phi[G.M - 1][j] - F.phi[G.M][j]) / G.hx;
        F.Ey[G.M][j] = 0.0;
    }

    for (int i = 0; i <= G.M; ++i)
    {
        F.Ey[i][0] = (F.phi[i][0] - F.phi[i][1]) / G.hy;
        F.Ex[i][0] = -P.d_phi; 

        F.Ey[i][G.N] = 0.0;
    }
    for (int i = 1; i < G.M; ++i)
    {
        F.Ex[i][G.N] = -d0x(F.phi, i, G.N, G.hx);
    }
}

void computeIy(Fields& F, const Grid& G, const Params& P)
{
    for (int i = 1; i < G.M; ++i)
    {
        for (int j = 1; j < G.N; ++j)
        {
            F.iy[i][j] = P.KI1 * F.C[i][j] * F.Ey[i][j] + P.KI2 * d0y(F.C, i, j, G.hy);
        }
    }

    for (int i = 0; i <= G.M; ++i)
    {
        F.iy[i][0] = P.KI1 * F.C[i][0] * F.Ey[i][0] + P.KI2 * (F.C[i][1] - F.C[i][0]) / G.hy;
        F.iy[i][G.N] = P.KI1 * F.C[i][G.N] * F.Ey[i][G.N];
    }

    for (int j = 1; j < G.N; ++j)
    {
        F.iy[0][j] = P.KI1 * F.C[0][j] * F.Ey[0][j] + P.KI2 * d0y(F.C, 0, j, G.hy);
        F.iy[G.M][j] = P.KI1 * F.C[G.M][j] * F.Ey[G.M][j] + P.KI2 * d0y(F.C, G.M, j, G.hy);
    }
}

void restoreIx(Fields& F, const Grid& G, const Params& P)
{
    for (int j = 1; j < G.N; ++j)
    {
        F.ix[0][j] = (P.KI1 * F.C[0][j] * F.Ex[0][j]) / (1.0 - P.KI2 * P.tA);
        for (int i = 1; i <= G.M; ++i)
        {
            F.ix[i][j] = F.ix[i - 1][j] - G.hx * d0y(F.iy, i, j, G.hy);
        }
    }

    for (int i = 0; i <= G.M; ++i)
        F.ix[i][0] = P.KI1 * F.C[i][0] * F.Ex[i][0];

    for (int i = 1; i < G.M; ++i)
        F.ix[i][G.N] = P.KI1 * F.C[i][G.N] * F.Ex[i][G.N] + P.KI2 * d0x(F.C, i, G.N, G.hx);
    F.ix[0][G.N] = P.KI1 * F.C[0][G.N] * F.Ex[0][G.N] / (1.0 - P.KI2 * P.tA);
    F.ix[G.M][G.N] = P.KI1 * F.C[G.M][G.N] * F.Ex[G.M][G.N] / (1.0 + P.KI2 * P.tK);
}

void solveT_ADI(Fields& F, const Grid& G, const Params& P)
{
    vector<vector<double>> Tstar = F.T;
    vector<vector<double>> Tnew = F.T;

    for (int j = 1; j < G.N; j++)
    {
        vector<double> A(G.M + 1, 0.0);
        vector<double> B(G.M + 1, 0.0);
        vector<double> D(G.M + 1, 0.0);
        vector<double> R(G.M + 1, 0.0);
        for (int i = 1; i < G.M; i++)
        {
            double Vx = F.Vx[i][j];
            double Vy = F.Vy[i][j];

            if (Vx >= 0.0)
            {
                A[i] = 1.0 / (P.Re * P.Pr * G.hx * G.hx) + Vx / G.hx;
                B[i] = -2.0 / G.ht - 2.0 / (P.Re * P.Pr * G.hx * G.hx) - Vx / G.hx;
                D[i] = 1.0 / (P.Re * P.Pr * G.hx * G.hx);
            }
            else
            {
                A[i] = 1.0 / (P.Re * P.Pr * G.hx * G.hx);
                B[i] = -2.0 / G.ht - 2.0 / (P.Re * P.Pr * G.hx * G.hx) + Vx / G.hx;
                D[i] = 1.0 / (P.Re * P.Pr * G.hx * G.hx) - Vx / G.hx;
            }
            double dCdy = upwindY(F.T, F.Vy, i, j, G.hy);
            double diffusionY = (F.T[i][j + 1] - 2.0 * F.T[i][j] + F.T[i][j - 1]) / (G.hy * G.hy);

            R[i] = -2.0 * F.T[i][j] / G.ht + Vy * dCdy - diffusionY / (P.Re * P.Pr) - P.BE * (F.ix[i][j] * F.Ex[i][j] + F.iy[i][j] * F.Ey[i][j]) / P.Pe;
        }

        A[0] = 0.0;
        B[0] = 1.0;
        D[0] = -1.0;
        R[0] = -G.hx * P.k * (F.T[0][j] - 1);

        A[G.M] = -1.0;
        B[G.M] = 1.0;
        D[G.M] = 0.0;
        R[G.M] = -G.hx * P.k * (F.T[G.M][j] - 1);
        vector<double> sol = thomasSolver(A, B, D, R);

        for (int i = 0; i <= G.M; i++)
            Tstar[i][j] = sol[i];
    }

    for (int i = 1; i < G.M; i++)
    {
        vector<double> A(G.N + 1, 0.0);
        vector<double> B(G.N + 1, 0.0);
        vector<double> D(G.N + 1, 0.0);
        vector<double> R(G.N + 1, 0.0);
        for (int j = 1; j < G.N; j++)
        {
            double Vx = F.Vx[i][j];
            double Vy = F.Vy[i][j];

            if (Vy >= 0.0)
            {
                A[j] = 1.0 / (P.Re * P.Pr * G.hy * G.hy) + Vy / G.hy;
                B[j] = -2.0 / G.ht - 2.0 / (P.Re * P.Pr * G.hy * G.hy) - Vy / G.hy;
                D[j] = 1.0 / (P.Re * P.Pr * G.hy * G.hy);
            }
            else
            {
                A[j] = 1.0 / (P.Re * P.Pr * G.hy * G.hy);
                B[j] = -2.0 / G.ht - 2.0 / (P.Re * P.Pr * G.hy * G.hy) + Vy / G.hy;
                D[j] = 1.0 / (P.Re * P.Pr * G.hy * G.hy) - Vy / G.hy;
            }
            double dCdx = upwindX(Tstar, F.Vx, i, j, G.hx);
            double diffusionX = (Tstar[i + 1][j] - 2.0 * Tstar[i][j] + Tstar[i - 1][j]) / (G.hx * G.hx);

            R[j] = -2.0 * Tstar[i][j] / G.ht + Vx * dCdx - diffusionX / (P.Re * P.Pr) - P.BE * (F.ix[i][j] * F.Ex[i][j] + F.iy[i][j] * F.Ey[i][j]) / P.Pe;
        }

        A[0] = 0.0;
        B[0] = 1.0;
        D[0] = 0.0;
        R[0] = 1.0;

        A[G.N] = -1.0;
        B[G.N] = 1.0;
        D[G.N] = 0.0;
        R[G.N] = 0.0;
        vector<double> sol = thomasSolver(A, B, D, R);

        for (int j = 0; j <= G.N; j++)
            Tnew[i][j] = sol[j];
    }
    F.T = Tnew;
    for (int j = 0; j <= G.N; ++j)
    {
        F.T[0][j] = F.T[1][j] - G.hx * P.k * (F.T[1][j] - 1.0);
        F.T[G.M][j] = F.T[G.M - 1][j] - G.hx * P.k * (F.T[G.M - 1][j] - 1.0);
    }
}

void solveXi_ADI(Fields& F, const Grid& G, const Params& P)
{
    vector<vector<double>> Xistar = F.xi;
    vector<vector<double>> Xinew = F.xi;

    for (int j = 1; j < G.N; j++)
    {
        vector<double> A(G.M + 1, 0.0);
        vector<double> B(G.M + 1, 0.0);
        vector<double> D(G.M + 1, 0.0);
        vector<double> R(G.M + 1, 0.0);
        for (int i = 1; i < G.M; i++)
        {
            double Vx = F.Vx[i][j];
            double Vy = F.Vy[i][j];

            if (Vx >= 0.0)
            {
                A[i] = 1.0 / (P.Re * G.hx * G.hx) + Vx / G.hx;
                B[i] = -2.0 / G.ht - 2.0 / (P.Re * G.hx * G.hx) - Vx / G.hx;
                D[i] = 1.0 / (P.Re * G.hx * G.hx);
            }
            else
            {
                A[i] = 1.0 / (P.Re * G.hx * G.hx);
                B[i] = -2.0 / G.ht - 2.0 / (P.Re * G.hx * G.hx) + Vx / G.hx;
                D[i] = 1.0 / (P.Re * G.hx * G.hx) - Vx / G.hx;
            }
            double dCdy = upwindY(F.xi, F.Vy, i, j, G.hy);
            double diffusionY = (F.xi[i][j + 1] - 2.0 * F.xi[i][j] + F.xi[i][j - 1]) / (G.hy * G.hy);
            double rhssin = P.Gr * ((P.AC * upwindY(F.C, F.Vy, i, j, G.hy) + P.AT * upwindY(F.T, F.Vy, i, j, G.hy)) * sin(P.gamma)) / pow(P.Re, 2);
            double rhscos = P.Gr * ((P.AC * upwindX(F.C, F.Vx, i, j, G.hx) + P.AT * upwindX(F.T, F.Vx, i, j, G.hx)) * cos(P.gamma)) / pow(P.Re, 2);
            R[i] = -2.0 * F.xi[i][j] / G.ht + Vy * dCdy - diffusionY / P.Re + rhssin - rhscos;
        }

        A[0] = 0.0;
        B[0] = 1.0;
        D[0] = 0.0;
        R[0] = 2 * (F.psi[1][j] - F.psi[0][j]) / pow(G.hx, 2);

        A[G.M] = 0.0;
        B[G.M] = 1.0;
        D[G.M] = 0.0;
        R[G.M] = 2 * (F.psi[G.M - 1][j] - F.psi[G.M][j]) / pow(G.hx, 2);
        vector<double> sol = thomasSolver(A, B, D, R);

        for (int i = 0; i <= G.M; i++)
            Xistar[i][j] = sol[i];
    }

    for (int i = 1; i < G.M; i++)
    {
        vector<double> A(G.N + 1, 0.0);
        vector<double> B(G.N + 1, 0.0);
        vector<double> D(G.N + 1, 0.0);
        vector<double> R(G.N + 1, 0.0);
        for (int j = 1; j < G.N; j++)
        {
            double Vx = F.Vx[i][j];
            double Vy = F.Vy[i][j];

            if (Vy >= 0.0)
            {
                A[j] = 1.0 / (P.Re * G.hy * G.hy) + Vy / G.hy;
                B[j] = -2.0 / G.ht - 2.0 / (P.Re * G.hy * G.hy) - Vy / G.hy;
                D[j] = 1.0 / (P.Re * G.hy * G.hy);
            }
            else
            {
                A[j] = 1.0 / (P.Re * G.hy * G.hy);
                B[j] = -2.0 / G.ht - 2.0 / (P.Re * G.hy * G.hy) + Vy / G.hy;
                D[j] = 1.0 / (P.Re * G.hy * G.hy) - Vy / G.hy;
            }
            double dCdx = upwindX(Xistar, F.Vx, i, j, G.hx);
            double diffusionX = (Xistar[i + 1][j] - 2.0 * Xistar[i][j] + Xistar[i - 1][j]) / (G.hx * G.hx);
            double rhssin = P.Gr * ((P.AC * upwindY(F.C, F.Vy, i, j, G.hy) + P.AT * upwindY(F.T, F.Vy, i, j, G.hy)) * sin(P.gamma)) / pow(P.Re, 2);
            double rhscos = P.Gr * ((P.AC * upwindX(F.C, F.Vx, i, j, G.hx) + P.AT * upwindX(F.T, F.Vx, i, j, G.hx)) * cos(P.gamma)) / pow(P.Re, 2);
            R[j] = -2.0 * Xistar[i][j] / G.ht + Vx * dCdx - diffusionX / P.Re + rhssin - rhscos;
        }

        A[0] = 0.0;
        B[0] = 1.0;
        D[0] = 0.0;
        R[0] = 0.0;
        //R[0] = 2 * (F.psi[i][1] - F.psi[i][0]) / pow(G.hy, 2);

        A[G.N] = 0.0;
        B[G.N] = 1.0;
        D[G.N] = 0.0;
        R[G.N] = 0.0;
        //R[G.N] = 2 * (F.psi[i][G.N - 1] - F.psi[i][G.N]) / pow(G.hy, 2);
        vector<double> sol = thomasSolver(A, B, D, R);

        for (int j = 0; j <= G.N; j++)
            Xinew[i][j] = sol[j];
    }
    F.xi = Xinew;
    for (int j = 0; j <= G.N; ++j)
    {
        F.xi[0][j] = 2 * (F.psi[1][j] - F.psi[0][j]) / pow(G.hx, 2);
        F.xi[G.M][j] = 2 * (F.psi[G.M - 1][j] - F.psi[G.M][j]) / pow(G.hx, 2);
    }
}

void solvePsi(Fields& F, const Grid& G, const Params& P)
{
    vector<vector<double>> psi_star = F.psi;
    vector<vector<double>> psi_new = F.psi;

    double err = 1.0;
    int iter = 0;
    int iterMax = 100;
    double K_psi = 2 * G.hx * G.hx / G.htau1;
    double beta = G.hx / G.hy;
    while (err > P.eps && iter < iterMax)
    {
        vector<vector<double>> psi_old = F.psi;
        for (int j = 1; j < G.N; ++j)
        {
            vector<double> A(G.M + 1), B(G.M + 1), D(G.M + 1), R(G.M + 1);
            for (int i = 1; i < G.M; ++i)
            {
                A[i] = 1.0;
                B[i] = -2 - K_psi;
                D[i] = 1.0;
                R[i] = -pow(beta, 2) * (F.psi[i][j - 1] - (2 - K_psi / pow(beta, 2)) * F.psi[i][j] + F.psi[i][j + 1]) + pow(G.hx, 2) * F.xi[i][j];
            }

            A[0] = 0; B[0] = 1; D[0] = 0; R[0] = 1;
            A[G.M] = 0; B[G.M] = 1; D[G.M] = 0; R[G.M] = 0;

            auto sol = thomasSolver(A, B, D, R);
            for (int i = 0; i <= G.M; ++i)
                psi_star[i][j] = sol[i];
        }

        for (int i = 1; i < G.M; ++i)
        {
            vector<double> A(G.N + 1), B(G.N + 1), D(G.N + 1), R(G.N + 1);

            for (int j = 1; j < G.N; ++j)
            {
                A[j] = 1.0;
                B[j] = -2 - K_psi / pow(beta, 2);
                D[j] = 1.0;
                R[j] = -(psi_star[i - 1][j] - (2 - K_psi) * psi_star[i][j] + psi_star[i + 1][j]) / pow(beta, 2) + pow(G.hy, 2) * F.xi[i][j];
            }

            A[0] = 0; B[0] = 1; D[0] = 0; R[0] = 1 - i * G.hx;
            A[G.N] = 0; B[G.N] = 1; D[G.N] = 0; R[G.N] = 1 - i * G.hx;

            auto sol = thomasSolver(A, B, D, R);

            for (int j = 0; j <= G.N; ++j)
                psi_new[i][j] = sol[j];
        }

        for (int j = 0; j <= G.N; ++j)
        {
            psi_new[0][j] = 1.0;
            psi_new[G.M][j] = 0;
        }

        err = 0.0;

        for (int i = 0; i <= G.M; ++i)
            for (int j = 0; j <= G.N; ++j)
                err = max(err, fabs(psi_new[i][j] - psi_old[i][j]));

        cout << "Iter = " << iter << " errPsi" << err << endl;
        F.psi = psi_new;
        iter++;
    }
}

void updateVelocity(Fields& F, const Grid& G, const Params& P)
{
    for (int i = 1; i < G.M; ++i)
    {
        for (int j = 1; j < G.N; ++j)
        {
            F.Vx[i][j] = d0y(F.psi, i, j, G.hy);
            F.Vy[i][j] = -d0x(F.psi, i, j, G.hx);
        }
    }

    for (int i = 0; i <= G.M; ++i)
    {
        F.Vx[i][0] = 0.0;
        F.Vx[i][G.N] = 0.0;

        F.Vy[i][0] = 1.0;
        F.Vy[i][G.N] = 1.0;
    }

    for (int j = 0; j <= G.N; ++j)
    {
        F.Vx[0][j] = 0.0;
        F.Vx[G.M][j] = 0.0;

        F.Vy[0][j] = 0.0;
        F.Vy[G.M][j] = 0.0;
    }
}