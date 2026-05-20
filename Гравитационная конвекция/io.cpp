#include "io.h"
#include "grid.h"
#include <fstream>

void saveCSV(const Fields& F, const std::string& filename) {
    std::ofstream file(filename);

    int M = F.C.size();
    int N = F.C[0].size();

    file << "i,j,C,T,phi,psi,Vx,Vy\n";

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            file << i << "," << j << ","
                << F.C[i][j] << ","
                << F.T[i][j] << ","
                << F.phi[i][j] << ","
                << F.psi[i][j] << ","
                << F.Vx[i][j] << ","
                << F.Vy[i][j] << "\n";
        }
    }

    file.close();
}

void saveField3D(const vector<vector<double>>& F,
    const Grid& G,
    const string& name)
{
    ofstream file(name);

    for (int i = 0; i <= G.M; ++i)
    {
        for (int j = 0; j <= G.N; ++j)
        {
            double x = i * G.hx;
            double y = j * G.hy;

            file << x << "," << y << "," << F[i][j] << "\n";
        }
    }
}