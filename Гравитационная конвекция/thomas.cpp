#include "thomas.h"

vector<double> thomasSolver(
    const vector<double>& A,
    const vector<double>& B,
    const vector<double>& D,
    const vector<double>& R
)
{
    int n = B.size();

    vector<double> c(n), d(n), x(n);

    // forward sweep
    c[0] = D[0] / B[0];
    d[0] = R[0] / B[0];

    for (int i = 1; i < n; i++) {

        double denom = B[i] - A[i] * c[i - 1];

        c[i] = (i < n - 1) ? D[i] / denom : 0.0;
        d[i] = (R[i] - A[i] * d[i - 1]) / denom;
    }

    // backward sweep
    x[n - 1] = d[n - 1];

    for (int i = n - 2; i >= 0; i--) {
        x[i] = d[i] - c[i] * x[i + 1];
    }

    return x;
}