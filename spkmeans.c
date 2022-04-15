
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

// need to modify this:
int main(int argc, char *argv[])
{
    int N, d, k, goal;
    FILE *input;
    double **output;

    if (!isValidInput(argc, argv, &input, &k, &goal))
    {
        printf("Invalid Input!\n");
        return 1;
    }

    data = read_matrix(input);

    output = execute_goal(data, n, d, *k, **mu, goal);

    return 0;
}

void read_matrix(FILE **input, double **matrix, double **mu, int n, int d, int k)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            if (j == d - 1)
            {
                fscanf(*input, "%le", &matrix[i][j]);
            }
            else
            {
                fscanf(*input, "%le,", &matrix[i][j]);
            }
            if (i < k)
            {
                mu[i][j] = matrix[i][j];
            }
        }
    }
}

void print_to_output(FILE **output, double **mu, int k, int d)
{
    int i, j;

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            fprintf(*output, "%.4f", mu[i][j]);
            if (j != d - 1)
            {
                fprintf(*output, ",");
            }
            else
            {
                fprintf(*output, "\n");
            }
        }
    }
    fclose(*output);
}

int isValidInput(int argc, char *argv[], FILE **input, FILE **output, int *k, int *max_iter)
{
    int no_max_iter_indent;
    /*chceck correct number of args*/
    if (argc != 5 && argc != 4)
    {
        return 0;
    }

    *k = isInt(argv[1]);

    if (argc == 5)
    {
        no_max_iter_indent = 1;
        *max_iter = isInt(argv[2]);
    }
    else
    {
        no_max_iter_indent = 0;
        *max_iter = 200; /*default*/
    }

    /*check first to args are ints*/
    if (!*k || !*max_iter)
    {
        return 0;
    }

    /*check if the files are "openable"*/
    *input = fopen(argv[2 + no_max_iter_indent], "r");
    *output = fopen(argv[3 + no_max_iter_indent], "w");

    if (!*input || !*output)
    {
        fclose(*input);
        fclose(*output);
        return 0;
    }
    return 1;
}

int isInt(char *intStr)
{
    int trueInt;
    int length = strlen(intStr);
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

int isDigit(char c)
{
    if ((c >= '0') && (c <= '9'))
        return 1;
    return 0;
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
