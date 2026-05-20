#include "operators.h"
using namespace std;

double d0x(
    const vector<vector<double>>& F,
    int i,
    int j,
    double hx
)
{
    return (F[i + 1][j] - F[i - 1][j]) / (2.0 * hx);
}
double d0y(
    const vector<vector<double>>& F,
    int i,
    int j,
    double hy
)
{
    return (F[i][j + 1] - F[i][j - 1]) / (2.0 * hy);
}
double upwindX(
    const vector<vector<double>>& F,
    const vector<vector<double>>& Vx,
    int i,
    int j,
    double hx
)
{
    if (Vx[i][j] >= 0.0)
        return (F[i][j] - F[i - 1][j]) / hx;
    else
        return (F[i + 1][j] - F[i][j]) / hx;
}
double upwindY(
    const vector<vector<double>>& F,
    const vector<vector<double>>& Vy,
    int i,
    int j,
    double hy
)
{
    if (Vy[i][j] >= 0.0)
        return (F[i][j] - F[i][j - 1]) / hy;
    else
        return (F[i][j + 1] - F[i][j]) / hy;
}
double laplacian(
    const vector<vector<double>>& F,
    int i,
    int j,
    double hx,
    double hy
)
{
    return
        (F[i + 1][j] - 2.0 * F[i][j] + F[i - 1][j]) / (hx * hx)
        + (F[i][j + 1] - 2.0 * F[i][j] + F[i][j - 1]) / (hy * hy);
}