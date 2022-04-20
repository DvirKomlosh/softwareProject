#include "spkmeans.h"
#include "algorithm.c"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define _POSIX_C_SOURCE 200809L

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
