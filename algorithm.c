#pragma once
#include "spkmeans.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define JAC_MAX_ITER 100
#define K_MAX_ITER 300
#define JAC_EPS 0.00001
#define K_EPS 0

#define sign(x) ((x >= 0) - (x < 0))

double *sort(double *array, int size)
{
    double temp;
    double *sorted = allocate_double_array(size);
    int i, j;

    for (i = 0; i < size; i++)
    {
        sorted[i] = array[i];
    }

    for (i = 0; i < size; ++i)
    {
        for (j = i + 1; j < size; ++j)
        {
            if (sorted[i] > sorted[j])
            {
                temp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = temp;
            }
        }
    }
    return sorted;
}

double **create_T(double **jacobi_mat, double *sorted_eigenvals, int k, int n)
{
    double **T;
    T = create_U(jacobi_mat, sorted_eigenvals, k, n);
    normalize(T, k, n);
    return T;
}

void normalize(double **U, int k, int n)
{
    int i, j;
    double norm;
    for (i = 0; i < n; i++)
    {
        norm = 0;
        for (j = 0; j < k; j++)
        {
            norm += U[i][j] * U[i][j];
        }
        norm = sqrt(norm);
        if (norm == 0) continue;
        for (j = 0; j < k; j++)
        {
            U[i][j] = U[i][j] / norm;
        }
    }
}

double **create_U(double **jacobi_mat, double *sorted_eigenvals, int k, int n)
{
    int i, j = 0, m;
    /* int time_to_stop_i = k - 1; */
    double current_eigenvalue;
    double **U;

    U = allocate_double_matrix(n, k);
    for (i = 0; i < k; i++)
    {
        current_eigenvalue = sorted_eigenvals[i];
        j = 0;
        while (current_eigenvalue != jacobi_mat[0][j])
        {
            j++;
        }
        /* if (i == time_to_stop_i) break; */
        for (m = 0; m < n ; m++)
        {
            /* the first row is the eigen values, so we skip it*/
            U[m][i] = jacobi_mat[m+1][j];
        }
    }
    return U;
}

double **allocate_double_matrix(int length, int width)
{
    int i;
    double **matrix = (double **)malloc(length * sizeof(double *));
    
    for (i = 0; i < length; i++)
    {
        matrix[i] = allocate_double_array(width);
    }
    return matrix;
}

double *allocate_double_array(int dim)
{
    return (double *)malloc(dim * sizeof(double));
}

/* gets the data points and n,d - the dim of the matrix,
   returns the weighted matrix */
double **wam(double **data, int n, int d)
{
    int i, j;

    /* allocating space for wam:*/
    double **wam_mat = allocate_double_matrix(n, n);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            wam_mat[i][j] = exp(-dist(data[i], data[j], d) / 2.0);
            if (i == j)
            {
                wam_mat[i][i] = 0;
            }
        }
    }

    return wam_mat;
}

/* gets the wam ,n dim of matrix
   returns the degrees diagonal matrix diagonal values*/
double *ddg(double **wam_mat, int n)
{
    int i, j;
    /* allocating space for ddg: */
    double *ddg_mat = allocate_double_array(n);

    for (i = 0; i < n; i++)
    {
        ddg_mat[i] = 0;
        for (j = 0; j < n; j++)
        {
            ddg_mat[i] += wam_mat[i][j];
        }
    }

    return ddg_mat;
}

/* gets the ddg,wam and n(the dim of matrix)
   returns lnorm matrix, destroys ddg in the process*/
double **lnorm(double **wam_mat, double *ddg_diag, int n)
{
    int i, j;
    double **lnorm_mat = allocate_double_matrix(n, n);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            lnorm_mat[i][j] = -wam_mat[i][j] * (1 / (sqrt(ddg_diag[i] * ddg_diag[j])));
            if (i == j)
                lnorm_mat[i][j] += 1;
        }
    }

    return lnorm_mat;
}

void get_rotation_values(int *i, int *j, double *c, double *s, double **A, int dim)
{
    int index_i, index_j;
    double max_off_diag = 0.0, t, theta;

    for (index_i = 0; index_i < dim; index_i++)
    {
        for (index_j = index_i + 1; index_j < dim; index_j++)
        {
            if (fabs(A[index_i][index_j]) > max_off_diag)
            {
                *i = index_i;
                *j = index_j;
                max_off_diag = fabs(A[index_i][index_j]);
            }
        }
    }

    theta = (A[*j][*j] - A[*i][*i]) / (2 * A[*i][*j]);
    t = (sign(theta)) / (fabs(theta) + sqrt(theta * theta + 1));
    *c = 1 / (sqrt(t * t + 1));
    *s = t * (*c);

}

/* changes A to A':
   used in jacobi */
void update_A(double **A, int i, int j, double c, double s, int dim)
{
    int index;
    double tempii, tempjj, ari, arj;
    for (index = 0; index < dim; index++)
    {
        if (index == i || index == j)
            continue;
        ari = A[index][i];
        arj = A[index][j];
        A[index][i] = c * ari - s * arj;
        A[i][index] = A[index][i];
        A[index][j] = c * arj + s * ari;
        A[j][index] = A[index][j];
    }
    tempii = A[i][i];
    tempjj = A[j][j];
    A[i][i] = c * c * tempii + s * s * tempjj - 2 * s * c * A[i][j];
    A[j][j] = s * s * tempii + c * c * tempjj + 2 * s * c * A[i][j];

    A[i][j] = 0;
    A[j][i] = 0;
}

/* updates V , used in jacobi */
void update_V(double **V, int i, int j, double c, double s, int dim)
{
    int index;
    double vri, vrj;

    for (index = 0; index < dim; index++)
    {
        vri = V[index][i];
        vrj = V[index][j];
        V[index][i] = c * vri - s * vrj;
        V[index][j] = s * vri + c * vrj;
    }
}

/* convergacne metric for Jacobi */
double off_diag_squared(double **matrix, int n)
{
    int i, j;
    double sum = 0;
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            sum += 2 * matrix[i][j] * matrix[i][j];
        }
    }
    return sum;
}

void set_to_identity(double **V, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            V[i][j] = (i == j);
        }
    }
}

/* returns the eigenvalues and eigenvectors of the matrix A,
   n dim of matrix, epsilon used to check convergence
   changes A in the process */
double **jacobi(double **A, int n, int max_iter, double epsilon)
{
    bool convarged = 0;
    double c, s, offA;
    int i, j, current_iter = 0;
    double **eigens;

    double **V = allocate_double_matrix(n, n);

    set_to_identity(V, n);
    offA = off_diag_squared(A, n);
    while (!convarged)
    {
        get_rotation_values(&i, &j, &c, &s, A, n);
        update_A(A, i, j, c, s, n);
        update_V(V, i, j, c, s, n);

        convarged = (fabs(off_diag_squared(A, n) - offA) < epsilon);
        offA = off_diag_squared(A, n);
        current_iter++;
        if (current_iter == max_iter)
        {
            convarged = true;
        }
    }
    eigens = allocate_double_matrix(n + 1, n);

    /* copies the eigenvalues and eigenvectors to the returned matrix */
    for (i = 0; i < n; i++)
    {
        eigens[0][i] = A[i][i];
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            eigens[i + 1][j] = V[i][j];
        }
    }

    free_matrix(V, n);
    return eigens;
}

double **diag_to_mat(double *diag, int n)
{
    int i, j;
    double **A;

    A = allocate_double_matrix(n, n);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = (i == j) * diag[i];
        }
    }
    return A;
}

/* gets sorted list of eigen values
   returns k acording to the eigengap method. */
int eigen_gap(double *eigen_values, int length)
{
    int i, max_index = 0;
    double max_eigen_gap = 0;
    if (length == 1) return 1;
    for (i = 0; i < (length / 2); i++)
    {
        if (fabs(eigen_values[i] - eigen_values[i + 1]) > max_eigen_gap)
        {
            max_eigen_gap = fabs(eigen_values[i] - eigen_values[i + 1]);
            max_index = i;
        }
    }
    return max_index + 1;
}

void Kmeans(double **matrix, double **mu, int n, int d, int k, int max_iter, double EPSILON)
{
    double minimal_distance;
    int curr_cluster, iter, j, i;
    int convergacne = 0;
    int *cluster_size = (int *)malloc(k * sizeof(int));
    double **mu_next = (double **)malloc(k * sizeof(double *));

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

        /* adding x to clusters */
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
    int i;
    double total = 0;
    for (i = 0; i < dim; i++)
    {
        total += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    return sqrt(total);
}

double **zeros(int n,int d)
{
    int j,i;
    double** T;
    T = allocate_double_matrix(n,d);
    for(i=0 ; i<n ; i++)
    {
        for(j=0 ; j<d ; j++)
        {
            T[i][j]=0.0;
        }
    }
    return T;
}

void print_to_output(double **output, int n, int d)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            printf("%.4f", output[i][j]);
            if (j != d - 1)
            {
                printf(",");
            }
            else
            {
                printf("\n");
            }
        }
    }
}

void free_matrix(double** matrix_to_free, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        free(matrix_to_free[i]);
    }
    free(matrix_to_free);
}

double **execute_goal(double **data, int n, int d, int *k, double **mu, int goal)
{
    double *sorted_eigenvals, *ddg_list_result;
    double **ddg_result, **wam_result, **lnorm_result, **jacobi_result, **T;

    if (goal == 5)
        return jacobi(data, n, JAC_MAX_ITER, JAC_EPS);

    if (goal == 1)
    {
        Kmeans(data, mu, n, d, *k, K_MAX_ITER, K_EPS);
        return mu;
    }

    wam_result = wam(data, n, d);
    if (goal == e_wam)
        return wam_result;

    ddg_list_result = ddg(wam_result, n);
    if (goal == e_ddg)
    {
        ddg_result = diag_to_mat(ddg_list_result, n);
        free(ddg_list_result);
        free_matrix(wam_result, n);
        return ddg_result;
    }

    lnorm_result = lnorm(wam_result, ddg_list_result, n);
    free_matrix(wam_result, n);
    free(ddg_list_result);
    if (goal == e_lnorm)
        return lnorm_result;

    if (goal == 6)
    {
        jacobi_result = jacobi(lnorm_result, n, JAC_MAX_ITER, JAC_EPS);
        sorted_eigenvals = sort(jacobi_result[0], n);
        if (*k == 0)
        {
            *k = eigen_gap(sorted_eigenvals, n);
        }
        T = create_T(jacobi_result, sorted_eigenvals, *k, n);
        free_matrix(lnorm_result, n);
        free_matrix(jacobi_result, n + 1);
        free(sorted_eigenvals);
        return T;
    }

    printf("Something very wrong happened.\n");
    return NULL;
}
