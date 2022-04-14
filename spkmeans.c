
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define sign(x) ((x >= 0) - (x < 0))

#define JAC_MAX_ITER 100
#define K_MAX_ITER 300
#define JAC_EPS 0.00001
#define K_EPS 0

int main(int argc, char *argv[]) // to do
{
}

//------------------------------------------------
// goal index:
// spk 1
// wam 2
// ddg 3
// lnorm 4
// jacobi 5
// kmeans 6
//------------------------------------------------

double **execute_goal(double **data, int n, int d, int *k, double **mu, int goal)
{
    double *sorted_eigenvals;

    if (goal == 5)
        return jacobi(data, n, JAC_MAX_ITER, JAC_EPS);

    if (goal == 1)
    {
        kmeans(data, mu, n, d, *k, K_MAX_ITER, K_EPS);
        return mu;
    }

    double **wam = wam(data, n, d);
    if (goal == 2)
        return wam;

    double *ddg = ddg(wam, n);
    if (goal == 3)
        ddg = diag_to_mat();
    return ddg;

    double **lnorm = lnorm(wam, ddg, n);
    if (goal == 4)
        return lnorm;

    if (goal == 6)
    {
        double **jacobi = jacobi(lnorm, n, JAC_MAX_ITER, JAC_EPS);
        sorted_eigenvals = sorted(jacobi[0], n);
        *k = eigen_gap(sorted_eigenvals, n);
        T = create_T(jacobi, sorted_eigenvals, *k);
        return T;
    }
}
