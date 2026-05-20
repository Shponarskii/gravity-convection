#include <iostream>
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "io.h"
#include "params.h"
#include "solvers.h"

using namespace std;

int main() {

    // размерные параметры
    double H_d = 0.007;
    double L_d = 0.02;
    double gamma = 0;
    double V0_d = 1e-4;
    double C0_d = 100;
    double T0_d = 293.15;
    double dphi_d = -0.4;
    double rho0_d = 1002.5;
    double nu_d = 1.006e-6;
    double D1_d = 1.33e-9;
    double D2_d = 2.05e-9;
    double alpha_d = 1.425e-7;
    double F_d = 96485.33838;
    double R_d = 8.31446;
    double cp_d = 4083;
    double AC_d = 0.04127;
    double AT_d = -0.23;
    double g_d = 9.80665;
    double n_1 = 1;
    double n_2 = 1;
    double z_1 = 1;
    double z_2 = -1;
    double D_d = D1_d * D2_d * (z_1 - z_2) / (z_1 * D1_d - z_2 * D2_d);
    double k_p = 0.6;
    double k_d = k_p / (alpha_d * rho0_d * cp_d);
    double delta_rho0_d = 4;
    double eps = 1e-3;
    double ht_d = 0.005;
    double t_max_d = 60;

    Params P;

    // безразмерные параметры
    double H = 1.0;
    double L = L_d / H_d;
    
    P.d_phi = (dphi_d * F_d) / (R_d * T0_d);

    P.Pe = (V0_d * H_d) / D_d;
    P.Re = (V0_d * H_d) / nu_d;   
    P.Pr = nu_d / alpha_d;
    P.Gr = (delta_rho0_d * g_d * pow(H_d, 3)) / (rho0_d * pow(nu_d, 2));

    P.D1 = D1_d / D_d;
    P.D2 = D2_d / D_d;
    P.KI1 = pow(z_1, 2) * n_1 * P.D1 + pow(z_2, 2) * n_2 * P.D2;
    P.KI2 = -(z_1 * n_1 * P.D1 + z_2 * n_2 * P.D2);

    P.Kphi = -P.KI2 / P.KI1;

    P.k = k_d * H_d;

    P.tA = 0.31;
    P.tK = 0.69;
    P.AC = C0_d * AC_d / delta_rho0_d;
    P.AT = T0_d * AT_d / delta_rho0_d;
    P.BE = (R_d * C0_d) / (rho0_d * cp_d);
    P.gamma = gamma;

    P.eps = eps;

    // --- параметры задачи ---
    int M = 100;
    int N = 100;
    double ht = ht_d * V0_d / H_d;
    double htau1 = ht;
    double htau2 = ht;

    // --- создание сетки ---
    Grid grid(M, N, L, ht, htau1, htau2);

    // --- поля ---
    Fields F(M, N);

    // --- начальные условия ---
    F.initialize();

    // --- time loop ---
    double t = 0.0;
    double t_max = t_max_d * V0_d / H_d;
    int step = 0;

    while (t < t_max)
    {
        cout << "======================================" << endl;

        cout << "Step: "
            << step
            << "   t = "
            << t << endl;

        // ==========================================
        // INNER ELECTRODIFFUSION ITERATIONS
        // ==========================================

        int couplingIter = 0;
        int couplingIterMax = 10;

        double residual = 1.0;

        while (residual > 1e-3 &&
            couplingIter < couplingIterMax)
        {
            auto Cold = F.C;
            auto PhiOld = F.phi;
            auto ExOld = F.Ex;
            auto EyOld = F.Ey;
            auto IxOld = F.ix;
            auto IyOld = F.iy;

            // ----------------------------------
            // ELECTRODIFFUSION BLOCK
            // ----------------------------------

            solveC_ADI(F, grid, P);
            double Cmin = 1e100;
            double Cmax = -1e100;

            for (int i = 0; i <= M; ++i)
            {
                for (int j = 0; j <= N; ++j)
                {
                    Cmin = min(Cmin, F.C[i][j]);
                    Cmax = max(Cmax, F.C[i][j]);
                }
            }

            cout << "Cmin = " << Cmin
                << "  Cmax = " << Cmax << endl;
            solvePhi_ADI(F, grid, P);

            computeE(F, grid, P);

            computeIy(F, grid, P);

            restoreIx(F, grid, P);

            // ----------------------------------
            // RESIDUALS ONLY FOR C AND PHI
            // ----------------------------------

            double errC = 0.0;
            double errPhi = 0.0;
            double errEx = 0.0;
            double errEy = 0.0;
            double errIx = 0.0;
            double errIy = 0.0;

            for (int i = 0; i <= M; ++i)
            {
                for (int j = 0; j <= N; ++j)
                {
                    errC = max(
                        errC,
                        fabs(F.C[i][j] - Cold[i][j])
                    );

                    errPhi = max(
                        errPhi,
                        fabs(F.phi[i][j] - PhiOld[i][j])
                    );
                    errEx = max(
                        errEx,
                        fabs(F.Ex[i][j] - ExOld[i][j])
                    );
                    errEy = max(
                        errEy,
                        fabs(F.Ey[i][j] - EyOld[i][j])
                    );
                    errIx = max(
                        errIx,
                        fabs(F.ix[i][j] - IxOld[i][j])
                    );
                    errIy = max(
                        errIy,
                        fabs(F.iy[i][j] - IyOld[i][j])
                    );
                }
            }

            residual = errC;

            couplingIter++;

            cout << "Electro iter = "
                << couplingIter
                << endl;

            cout << "errC   = "
                << errC
                << endl;

            cout << "errPhi = "
                << errPhi
                << endl;

            cout << "errEx = "
                << errEx
                << endl;

            cout << "errEy = "
                << errEy
                << endl;

            cout << "errIx = "
                << errIx
                << endl;

            cout << "errIy = "
                << errIy
                << endl;

            cout << "----------------------------------"
                << endl;
        }

        // ==========================================
        // HYDRODYNAMIC UPDATE
        // ==========================================
        auto TOld = F.T;
        auto XiOld = F.xi;
        auto PsiOld = F.psi;
        auto VxOld = F.Vx;
        auto VyOld = F.Vy;
        solveT_ADI(F, grid, P);

        solveXi_ADI(F, grid, P);

        solvePsi(F, grid, P);

        updateVelocity(F, grid, P);

        double errT = 0.0;
        double errXi = 0.0;
        double errPsi = 0.0;
        double errVx = 0.0;
        double errVy = 0.0;

        double maxV = 0.0;

        for (int i = 0; i <= M; ++i)
        {
            for (int j = 0; j <= N; ++j)
            {
                errT = max(
                    errT,
                    fabs(F.T[i][j] - TOld[i][j])
                );

                errXi = max(
                    errXi,
                    fabs(F.xi[i][j] - XiOld[i][j])
                );

                errPsi = max(
                    errPsi,
                    fabs(F.psi[i][j] - PsiOld[i][j])
                );

                errVx = max(
                    errVx,
                    fabs(F.Vx[i][j] - VxOld[i][j])
                );

                errVy = max(
                    errVy,
                    fabs(F.Vy[i][j] - VyOld[i][j])
                );

                maxV = max(maxV, fabs(F.Vx[i][j]));
                maxV = max(maxV, fabs(F.Vy[i][j]));
            }
        }
        cout << "errT   = "
            << errT
            << endl;

        cout << "errXi  = "
            << errXi
            << endl;

        cout << "errPsi = "
            << errPsi
            << endl;

        cout << "errVx  = "
            << errVx
            << endl;

        cout << "errVy  = "
            << errVy
            << endl;

        cout << "maxV   = "
            << maxV
            << endl;
        // ==========================================
        // NEXT TIME STEP
        // ==========================================

        t += ht;

        step++;
    }

     //--- вывод ---
    saveField3D(F.C, grid, "C.csv");
    saveField3D(F.T, grid, "T.csv");
    saveField3D(F.phi, grid, "Phi.csv");
    saveField3D(F.Vx, grid, "Vx.csv");
    saveField3D(F.Vy, grid, "Vy.csv");
    saveField3D(F.Ex, grid, "Ex.csv");
    saveField3D(F.Ey, grid, "Ey.csv");
    saveField3D(F.psi, grid, "Psi.csv");
    saveField3D(F.xi, grid, "Xi.csv");

    std::cout << "Finished." << std::endl;

    return 0;
}