#pragma once
#include <stdio.h>

/* The requested goal enum */
enum goal_enum
{
    e_spk = 1,
    e_wam = 2,
    e_ddg = 3,
    e_lnorm = 4,
    e_jacobi = 5,
    e_kmeans = 6
};

/* Functions from spkmeans.c */
void read_matrix(FILE **input, double **matrix, int n, int d);
int isValidInput(int argc, char *argv[], FILE **input, enum goal_enum *goal);
void get_sizes(FILE **input, int *n, int *d);
int check_symmetry(double **data, int n, int d);
void print_to_output(double **output, int n, int d);

/* Functions from algorithm.c */

/* Allocation functions */
double **allocate_double_matrix(int length, int width);
double *allocate_double_array(int dim);
void free_matrix(double **matrix_to_free, int n);

/* Algorithmical functions */
double **wam(double **data, int n, int d);
double *ddg(double **wam_mat, int n);
double **lnorm(double **wam_mat, double *ddg_mat, int n);
double **jacobi(double **A, int n, int max_iter, double epsilon);
void Kmeans(double **matrix, double **mu, int n, int d, int k, int max_iter, double EPSILON);
double **execute_goal(double **data, int n, int d, int *k, double **mu, int goal);
double **diag_to_mat(double *diag, int n);

/* Inner functions */
void normalize(double **U, int k, int n);
double **create_U(double **jacobi_mat, double *sorted_eigenvals, int k, int n);
int isDigit(char c);
void set_to_identity(double **V, int n);
double dist(double *x1, double *x2, int dim);