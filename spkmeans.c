#include "spkmeans.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define _POSIX_C_SOURCE 200809L
#define JAC_MAX_ITER 100
#define K_MAX_ITER 300
#define JAC_EPS 0.00001
#define K_EPS 0

/* main function of the C side,
can compute the goals jacobi, wam, ddg,lnorm */
int main(int argc, char *argv[])
{
    int n, d;
    int k = -1; /* Placeholder variable, not used from the c main.*/
    enum goal_enum goal;
    FILE *input;
    double **data, **output;
    double **mu = NULL; /* placeholder pointer,not used from the c main*/

    /* check invalid input */
    if (!isValidInput(argc, argv, &input, &goal))
        return 1;

    get_sizes(&input, &n, &d);

    data = allocate_double_matrix(n, d);
    read_matrix(&input, data, n, d);

    if (goal == e_jacobi && !check_symmetry(data, n, d))
    {
        /* Isymmetric matrix given for Jacobi */
        printf("Invalid Input!\n");
        return 1;
    }
    /* out put should have the output matrix after executing the goal*/
    output = execute_goal(data, n, d, &k, mu, goal);

    if (goal == e_jacobi)
    {
        print_to_output(output, n + 1, n);
        free_matrix(output, n + 1);
    }
    else
    {
        print_to_output(output, n, n);
        free_matrix(output, n);
    }
    free_matrix(data, n);
    return 0;
}

/*  reads a file that holds a matrix of floats- of some
    places the values inside matrix                    */
void read_matrix(FILE **input, double **matrix, int n, int d)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            if (j == d - 1)
            {
                if (!fscanf(*input, "%le", &matrix[i][j]))
                {
                    /*file is not in the correct format*/
                    printf("An Error Has Occurred\n");
                    return;
                }
            }
            else
            {
                if (!fscanf(*input, "%le,", &matrix[i][j]))
                {
                    printf("An Error Has Occurred\n");
                    return;
                }
            }
        }
    }
}

/*  Validates the input, printing "Invalid Input!" if invalid
    returns 1 if valid, and parses the arguments while doing so
    (placing the FILE into input, goal into goal)              */
int isValidInput(int argc, char *argv[], FILE **input, enum goal_enum *goal)
{
    char *goal_string, *filename, *suffix;
    if (argc != 3)
    {
        /* Invalid amount of parameters */
        printf("Invalid Input!\n");
        return 0;
    }

    goal_string = argv[1];
    if (strcmp(goal_string, "wam") == 0)
        *goal = e_wam;
    else if (strcmp(goal_string, "ddg") == 0)
        *goal = e_ddg;
    else if (strcmp(goal_string, "lnorm") == 0)
        *goal = e_lnorm;
    else if (strcmp(goal_string, "jacobi") == 0)
        *goal = e_jacobi;
    else
    {
        /* Invalid value of the first parameter */
        printf("Invalid Input!\n");
        return 0;
    }

    /* check for thr right suffix of output */
    filename = argv[2];
    suffix = strrchr(filename, '.');
    if (suffix && (!strcmp(suffix, ".csv") || !strcmp(suffix, ".txt")))
        *input = fopen(filename, "r");
    if (*input == NULL)
    {
        /* Unopenable file */
        printf("Invalid Input!\n");
        return 0;
    }
    return 1;
}

/* computes the dimentions of the data matrix N,d
   from the input file, and assign them.
   printing an error message in case of invalid input*/
void get_sizes(FILE **input, int *N, int *d)
{
    int dtemp, Ntemp, found_d;
    int c;
    c = 0;
    found_d = 0;
    dtemp = 1;
    Ntemp = 0;
    if (feof(*input))
    {
        printf("An Error Has Occurred\n");
        return;
    }
    while (1)
    {
        c = fgetc(*input);

        if (c == '\n')
        {
            Ntemp++;
            if (!found_d)
            {
                *d = dtemp;
                found_d = 1;
            }
        }
        if (c == ',')
        {
            dtemp++;
        }

        if (feof(*input))
        {
            break;
        }
    }
    *N = Ntemp;
    rewind(*input);

    return;
}

/* checks if the matrix data is symetrical,
   returns 1 if it is, 0 otherwise         */
int check_symmetry(double **data, int n, int d)
{
    int i, j;
    if (n != d)
        return 0;
    for (i = 0; i < n; i++)
    {
        for (j = i; j < d; j++)
        {
            if (data[i][j] != data[j][i])
                return 0;
        }
    }
    return 1;
}

/*  Checks if a char array represents an int,
    returns the int, 0 if it does not        */
int isInt(char *intStr)
{
    int trueInt;
    int length = (int)strlen(intStr);
    int i;
    for (i = 0; i < length; i++)
    {
        if (!isDigit(intStr[i]))
            return 0;
    }
    trueInt = (int)strtol(intStr, (char **)NULL, 10);

    if (trueInt < 1)
    {
        return 0;
    }
    return trueInt;
}

/*  Checks if a char is a digit (0-9),
    true if it is                       */
int isDigit(char c)
{
    if ((c >= '0') && (c <= '9'))
        return 1;
    return 0;
}

/* a macro to compute the sign of a double,
   (we define sign(0) = 1)                   */
#define sign(x) ((x >= 0) - (x < 0))

/*  sortes a copy of array of size n, returns a pointer to the array*/
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

/*  creates the matrix T from the eigen values
    and the corresponding eigenvectors          */
double **create_T(double **eigen_mat, double *sorted_eigenvals, int k, int n)
{
    double **T;

    T = create_U(eigen_mat, sorted_eigenvals, k, n); /* creating the U matrix*/
    normalize(T, k, n);                              /* normaizing U to get T*/
    return T;
}

/* normalizes the rows of the matrix U, with dimentions n*d */
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
        /*in the case that the norm is 0, we do not change the vector,
          so we do not want to divide by 0*/
        if (norm == 0)
            continue;
        for (j = 0; j < k; j++)
        {
            U[i][j] = U[i][j] / norm;
        }
    }
}
/*  creates the U matrix, consisting of the first k
    eigenvectors of the lnorm matrix, sorted by eigen value */
double **create_U(double **jacobi_mat, double *sorted_eigenvals, int k, int n)
{
    int i, j = 0, m;
    double current_eigenvalue;
    double **U;

    U = allocate_double_matrix(n, k);
    for (i = 0; i < k; i++)
    {
        current_eigenvalue = sorted_eigenvals[i];
        j = 0;
        /*find the column index of the eigenvector
          corresponding to the current eigen value method*/
        while (current_eigenvalue != jacobi_mat[0][j])
        {
            j++;
        }
        for (m = 0; m < n; m++)
        {
            /* the first row is the eigen values, so we skip it*/
            U[m][i] = jacobi_mat[m + 1][j];
        }
    }
    return U;
}

/*alocates the space of a length * width space in the memory as a matrix*/
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

/*alocates an array of size dim in the memory*/
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

/* gets the ddg,wam and n (the dim of matrix)
   returns lnorm matrix, destroys ddg in the process*/
double **lnorm(double **wam_mat, double *ddg_diag, int n)
{
    int i, j;
    double **lnorm_mat = allocate_double_matrix(n, n);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            lnorm_mat[i][j] = -wam_mat[i][j] *
                              (1 / (sqrt(ddg_diag[i] * ddg_diag[j])));
            if (i == j)
                lnorm_mat[i][j] += 1;
        }
    }

    return lnorm_mat;
}

/* a helper method for Jacobi, placing the
   correct i,j,c,s into their adresses     */
void get_rotation_values(int *i, int *j, double *c,
                         double *s, double **A, int dim)
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

/* a helper method for Jacobi
   changes A to A' (one iteration of the jacobi algorithm)*/
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

/* a helper method for Jacobi
   updates V (one iteration of the jacobi algorithm)*/
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

/*  convergacne metric for Jacobi, returns the sum of squared
    elements on the off-diagonal of the square matrix of dimentions n*n */
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

/*  set an n*n square matrix to an identity matrix*/
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

        /*check convergence:*/
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

/* creates a diagonal n*n matrix from the array diag (of length n)
   with the diag[i]=matrix[i][i]                                   */
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
    if (length == 1)
        return 1;
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

/*preforms the k means algorithm on the points in matrix
  with the initilization points given inside mu         */
void Kmeans(double **matrix, double **mu, int n, int d,
            int k, int max_iter, double EPSILON)
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

/* computes the eulerian distance between vectors x1,x2
    of dimention dim                                    */
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

/*prints a matrix of size n*d to output
  with the conventions given for the project */
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

/* frees a matrix built of n rows */
void free_matrix(double **matrix_to_free, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        free(matrix_to_free[i]);
    }
    free(matrix_to_free);
}

/*  main function, executes the goal called from either python or C
    returns the matrix to be outputed for each goal
    (in the case goal = e_kmeans- an intermidate step to return to python
     before performing the full kmeans algorithm, we return the T matrix) */
double **execute_goal(double **data, int n, int d,
                      int *k, double **mu, int goal)
{
    double *sorted_eigenvals, *ddg_list_result;
    double **ddg_result, **wam_result, **lnorm_result, **jacobi_result, **T;

    if (goal == e_jacobi)
        return jacobi(data, n, JAC_MAX_ITER, JAC_EPS);

    if (goal == e_spk)
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

    if (goal == e_kmeans)
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
