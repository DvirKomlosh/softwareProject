#include "spkmeans.h"
#include "algorithm.c"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define _POSIX_C_SOURCE 200809L

int main(int argc, char *argv[])
{
    int n, d;
    int k = -1; /* Placeholder variable, as it's not used in a goal other than spk and kmeans. */
    enum goal_enum goal;
    FILE *input;
    double **data, **output;
    double **mu = NULL; /* placeholder pointer, we dont use it from the c main*/

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
                if (!fscanf(*input, "%le", &matrix[i][j]))
                {
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
                return 0;
        }
    }
    return 1;
}

int isInt(char *intStr)
{
    int trueInt;
    int length = (int) strlen(intStr);
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
