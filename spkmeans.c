
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define sign(x) ((x >= 0) - (x < 0))

int main(int argc, char *argv[]) // to do
{
}

// gets the data points and n,d - the dim of the matrix,
// returns the weighted matrix
double **wam(double **data, int n, int d)
{
    int i, j;

    // allocating space for wam:
    double **wam = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
    {
        wam[i] = (double *)malloc(n * sizeof(double));
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            wam[i][j] = exp(-dist(data[i], data[j], d) / 2.0);
        }
    }

    return wam;
}

// gets the wam ,n dim of matrix
// returns the degrees diagonal matrix diagonal values
double *ddg(double **wam, int n)
{
    // allocating space for ddg:
    double *ddg = (double *)malloc(n * sizeof(double));

    for (i = 0; i < n; i++)
    {
        ddg[i] = 0;
        for (j = 0; j < n; j++)
        {
            ddg[i] += wam[i][j];
        }
    }

    return ddg;
}

// void matrix_vector_multiplication(double *vector, double **matrix, double *res, int dim)
// {
//     int i, j;

//     for (i = 0; i < dim; i++)
//     {
//         res[i] = 0;
//         for (j = 0; j < dim; j++)
//         {
//             res[i] += matrix[i][j] * vector[j];
//         }
//     }
// }

// void matrix_multiplication(double **mat1, double **mat2, double **res, int dim)
// {
//     int i;
//     double **res = (double **)malloc(dim * sizeof(*double));
//     for (i = 0; i < dim; i++)
//     {
//         matrix_vector_multiplication(mat1[i], mat2, res[i], dim);
//     }
// }

// void one_over_sqrt_matrix(double **matrix, int dim)
// {
//     int i;
//     double **res = (double **)malloc(n * sizeof(double *));
//     for (i = 0; i < n; i++)
//     {
//         res[i] = (double *)malloc(n * sizeof(double));
//     }

//     for (i = 0; i < dim; i++)
//     {
//         res[i][i] = 1 / sqrt(res[i][i]);
//     }

//     return res;
// }

// gets the ddg,wam and n(the dim of matrix)
// returns lnorm matrix, destroys ddg in the process
double **lnorm(double **wam, double *ddg, int n)
{

    double **lnorm = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
    {
        lnorm[i] = (double *)malloc(n * sizeof(double));
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            lnorm[i][j] = -wam[i][j] * (1 / (sqrt(ddg[i] * ddg[j])));
            if (i == j)
                lnorm[i][j] += 1;
        }
    }

    return lnorm;
}

get_rotation_values(int *i, int *j, double *c, double *s, double **A, int dim)
{
    int index_i, index_j;
    double max_off_diag = 0.0, t, theta, c_temp, s_temp;

    for (index_i = 0; index_i < dim; index_i++)
    {
        for (index_j = index_i + 1; index_j < dim; index_j++)
        {
            if (abs(A[index_i][index_j]) > max_off_diag)
            {
                *i = index_i;
                *j = index_j;
                max_off_diag = abs(A[index_i][index_j]);
            }
        }
    }

    theta = (A[*j][*j] - A[*i][*i]) / (2 * A[*i][*j]);
    t = (sign(theta)) / (abs(theta) + sqrt(theta * theta + 1));
    c_temp = 1 / (sqrt(t * t + 1));
    s_temp = t * c;
    *c = c_temp;
    *s = s_temp;
}

// changes A to A':
void update_A(double **A, int i, int j, double c, double s, int dim)
{
    int index;

    for (index = 0; index_i < dim; index_i++)
    {
        A[index][i] = c * A[index][i] - s * A[index][j];
        A[i][index] = [index][i];
        A[index][j] = c * A[index][i] + s * A[index][j];
        A[j][index] = A[index][j];
    }
    A[i][i] = c * c * A[i][i] + s * s * A[j][j] - 2 * s * c * A[i][j];
    A[j][j] = s * s * A[i][i] + c * c * A[j][j] + 2 * s * c * A[i][j];
}

void update_V(double **V, int i, int j, double c, double s, int dim)
{
    int index;

    for (index = 0; index_i < dim; index_i++)
    {
        V[index][i] = c * V[index][i] - s * A[index][j];
        V[index][j] = c * V[index][i] + s * A[index][j];
    }
}

double off_diag_squared(double **matrix, n)
{
    int i, j;
    double sum;
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            sum += 2 * A[i][j] * A[i][j];
        }
    }
    return sum;
}

// returns the eigenvalues and eigenvectors of the matrix
// changes A in the process
double **jacobi(double **A, int n, double epsilon)
{
    bool convarged = false;
    double c, s, offA;
    int i, j;

    double **V = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
    {
        V[i] = (double *)malloc(n * sizeof(double));
    }

    set_to_identity(V);
    offA = off_diag_squared(A);
    while (!convarged)
    {

        get_rotation_values(A, &i, &j, &c, &s, n);
        update_A(A, i, j, c, s, n);
        update_V(V, i, j, c, s, n);

        convarged = ((off_diag_squared(A) - offA) < epsilon)
    }

    // A's diag should have all of the eigenvalues, V's rows should be the correspondant vectors.
    // couldn't understand how to output.
}

// gets list of eigen values
// returns k acording to the eigengap method.
int eigen_gap(double *eigen_values, int length)
{
    int i, max_index;
    double max_eigen_gap = 0;
    for (i = 1; i <= (length / 2); i++)
    {
        if (abs(eigen_values[i] - eigen_values[i + 1]) > max_eigen_gap)
        {
            max_eigen_gap = abs(eigen_values[i] - eigen_values[i + 1]);
            max_index = i;
        }
    }
    return max_index;
}

void Kmeans(double **matrix, double **mu, int n, int d, int k, int max_iter, double EPSILON)
{
    double minimal_distance;
    int curr_cluster, iter, j, i;
    int convergacne = 0;
    int *cluster_size = (int *)malloc(k * sizeof(int));
    double *mu_next = (double *)malloc(k * sizeof(double *));

    for (i = 0; i < k; i++)
    {
        mu_next[i] = (double *)malloc(d * sizeof(double));
    }

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            mu_next[i][j] = 0.0;
        }
    }

    /* Main Algorithem:  */
    for (iter = 0; !convergacne && iter < max_iter; iter++)
    {
        for (i = 0; i < k; i++)
        {
            cluster_size[i] = 0;
        }

        / adding x to clusters /
            for (i = 0; i < n; i++)
        {
            minimal_distance = dist(matrix[i], mu[0], d);
            curr_cluster = 0;
            for (j = 1; j < k; j++)
            {
                if (dist(matrix[i], mu[j], d) < minimal_distance)
                {
                    minimal_distance = dist(matrix[i], mu[j], d);
                    curr_cluster = j;
                }
            }

            cluster_size[curr_cluster] = cluster_size[curr_cluster] + 1;
            for (j = 0; j < d; j++)
            {
                mu_next[curr_cluster][j] += matrix[i][j];
            }
        }

        for (i = 0; i < k; i++)
        {
            for (j = 0; j < d; j++)
            {
                if (cluster_size[i] != 0)
                {
                    mu_next[i][j] = mu_next[i][j] / (double)cluster_size[i];
                }
            }
        }
        convergacne = 1;
        for (j = 0; j < k; j++)
        {
            if (dist(mu_next[j], mu[j], d) >= EPSILON)
            {
                convergacne = 0;
                break;
            }
        }
        for (i = 0; i < k; i++)
        {
            for (j = 0; j < d; j++)
            {
                mu[i][j] = mu_next[i][j];
                mu_next[i][j] = 0.0;
            }
        }
    }

    for (i = 0; i < k; i++)
    {
        free(mu_next[i]);
    }
    free(mu_next);
    free(cluster_size);
}

double dist(double *x1, double *x2, int dim)
{
    return sqrt(dist_squared);
}

double dist_squared(double *x1, double *x2, int dim)
{
    int i;
    double total = 0;
    for (i = 0; i < dim; i++)
    {
        total += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    return sqrt(total);
}