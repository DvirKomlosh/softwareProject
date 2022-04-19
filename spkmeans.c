#include "spkmeans.h"
#include "algorithm.c"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define _POSIX_C_SOURCE 200809L

#define JAC_MAX_ITER 100
#define K_MAX_ITER 300
#define JAC_EPS 0.00001
#define K_EPS 0

int main(int argc, char *argv[])
{
    int n, d;
    int k = -1; /* Placeholder variable, as it's not used in a goal other than spk and kmeans. */
    enum goal_enum goal;
    FILE *input;
    double **data, **output;
    double **mu = NULL; /* placeholder pointer, we dont use it from the c*/

    /* check invalid input */
    if (!isValidInput(argc, argv, &input, &goal))
    {
        printf("Invalid Input!\n");
        return 1;
    }
    get_sizes(&input, &n, &d);

    data = allocate_double_matrix(n, d);
    read_matrix(&input, data, n, d);

    if (goal == e_jacobi && !check_symmetry(data, n, d))
    {
        /* Isymmetric matrix given for Jacobi */
        printf("Invalid Input!\n");
        return 1;
    }
    output = execute_goal(data, n, d, &k, mu, goal);

    if (goal == e_jacobi)
    {
        print_to_output(output, n + 1, n);
    }
    else
        print_to_output(output, n, n);

    return 0;
}

void read_matrix(FILE **input, double **matrix, int n, int d)
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
        }
    }
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
/*
void get_sizes(FILE **input, int *n, int *d)
{
    char *lineptr;
    int line_read = 1, i;
    n = 0;

    while (line_read)
    {
        printf("read line");
        lineptr = NULL;
        if (fscanf(*input, "%[^\n]", lineptr) == 0)
        {
            line_read = 0;
            *d = 1;
            for (i = 0; i < (int)strlen(lineptr); i++)
            {
                if (lineptr[i] == ',')
                    d += 1;
            }
        }
        n += 1;
    }
    fseek(*input, 0, SEEK_SET);
    printf("done get_sizes\n");
}
*/

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
        printf("file was empty\n");
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
                return 1;
        }
    }
    return 0;
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

/*------------------------------------------------
 goal index:
 spk 1
 wam 2
 ddg 3
 lnorm 4
 jacobi 5
 kmeans 6 (an intermidate step, that returns T to the python)
--------------------------------------------------*/

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
        return ddg_result;
    }

    lnorm_result = lnorm(wam_result, ddg_list_result, n);
    if (goal == e_lnorm)
        return lnorm_result;

    if (goal == 6)
    {
        jacobi_result = jacobi(lnorm_result, n, JAC_MAX_ITER, JAC_EPS);
        sorted_eigenvals = sort(jacobi_result[0], n);
        *k = eigen_gap(sorted_eigenvals, n);
        T = create_T(jacobi_result, sorted_eigenvals, *k, n);
        return T;
    }

    printf("Something very wrong happened.\n");
    return NULL;
}
