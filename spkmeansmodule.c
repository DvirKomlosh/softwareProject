#include "spkmeans.h"
#include "algorithm.c"

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

PyMODINIT_FUNC PyInit_spkmeans(void);
static PyObject* fit(PyObject *self, PyObject *args);

/* This part of the code is taken mostly from the presentations of the course */
static PyMethodDef kmeansMethods[] = {
    {"fit", 
    (PyCFunction) fit,
    METH_VARARGS, 
    PyDoc_STR("The main algorithmical function, which receives a goal flag and performs the wanted goal")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduleDef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans",
    NULL, 
    -1, 
    kmeansMethods
};

PyMODINIT_FUNC PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduleDef);
    if (!m)
    {
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }
    return m;
}
/* End of current part */


/*
 * This function is the main function of this code file.
 * It receives all important parameters from the Python program.
 * This includes a goal flag, which represents the chosen goal.
 * Then, it activates the wanted function from the C file, 
 * and returns its result back to the Python program.
 * Information about PyObject and List Objects 
 * is taken from the Python documentations and Stack Overflow.
 * Most useful page - https://docs.python.org/3/c-api/list.html.
 */
static PyObject* fit(PyObject *self, PyObject *args)
{
    /* Declerations of useful variables */
    int initialize_i;
    int i, j;
    enum goal_enum goal;
    PyObject* po_primary_i;
    PyObject* po_mu_arr_i;
    PyObject* po_primary_i_j;
    PyObject* po_mu_arr_i_j;
    PyObject* po_return_mat_i;

    int n, d, k, goal_num;
    PyObject* po_primary;
    double** primary;
    PyObject* po_mu_arr;
    double** mu_arr = NULL;
    PyObject* po_return_mat;
    double** returned_mat;

    /* Receive the useful information from the user */
    printf("0\n");
    if (!PyArg_ParseTuple(args, "OiiiOi", &po_primary, &n, &d, &k, 
    &po_mu_arr, &goal_num))
    {
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }

    printf("0.5\n");
    if (goal_num == 1) goal = e_spk;
    else if (goal_num == 2) goal = e_wam;
    else if (goal_num == 3) goal = e_ddg;
    else if (goal_num == 4) goal = e_lnorm;
    else if (goal_num == 5) goal = e_jacobi;
    else if (goal_num == 6) goal = e_kmeans;
    else
    {
        /* Unexpected goal value */
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }

    /* Transform the given PyObjects to double matrices */
    printf("1\n");
    primary = (double **)malloc(n * sizeof(double *));
    if (!primary)
    {
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }
    printf("2\n");
    for (initialize_i = 0; initialize_i < n; initialize_i++)
    {
        primary[initialize_i] = (double *)malloc(d * sizeof(double));
        if (!primary[initialize_i])
        {
            printf("An Error Has Occurred\n");
            return Py_BuildValue("");
        }
    }
    
    printf("3\n");
    for (i = 0; i < n; i++)
    {
        po_primary_i = PyList_GetItem(po_primary, i);
        for (j = 0; j < d; j++)
        {
            po_primary_i_j = PyList_GetItem(po_primary_i, j);
            primary[i][j] = (double)PyFloat_AsDouble(po_primary_i_j);
        }
    }
    printf("4\n");
    if (goal == e_spk)
    {
        mu_arr = (double **)malloc(k * sizeof(double *));
        printf("5\n");
        if (!mu_arr)
        {
            printf("An Error Has Occurred\n");
            return Py_BuildValue("");
        }
        printf("6\n");
        for (initialize_i = 0; initialize_i < k; initialize_i++)
        {
            mu_arr[initialize_i] = (double *)malloc(d * sizeof(double));
            if (!mu_arr[initialize_i])
            {
                printf("An Error Has Occurred\n");
                return Py_BuildValue("");
            }
        }
        printf("7\n");
        for (i = 0; i < k; i++)
        {
            po_mu_arr_i = PyList_GetItem(po_mu_arr, i);
            for (j = 0; j < d; j++)
            {
                po_mu_arr_i_j = PyList_GetItem(po_mu_arr_i, j);
                mu_arr[i][j] = (double)PyFloat_AsDouble(po_mu_arr_i_j);
            }
        }
    }

    /* Activate the main C function */
    returned_mat = execute_goal(primary, n, d, &k, mu_arr, goal);
    printf("n1: %d\n", n);
    printf("d1: %d\n", d);
    printf("k1: %d\n", k);
    /* Transform the returned matrix to a PyObject (PyList) */   
    if (goal == e_spk) po_return_mat = PyList_New(k);
    else if (goal == e_jacobi) po_return_mat = PyList_New(n+1);
    else po_return_mat = PyList_New(n);
    printf("9\n");
    for (i = 0; i < PyList_Size(po_return_mat); i++)
    {
        if (goal == e_spk) po_return_mat_i = PyList_New(d);
        else if (goal == e_kmeans) po_return_mat_i = PyList_New(k);
        else po_return_mat_i = PyList_New(n);
        for (j = 0; j < PyList_Size(po_return_mat_i); j++)
        {
            PyList_SetItem(po_return_mat_i, j,
                Py_BuildValue("d", returned_mat[i][j]));
            printf("i = %d, j = %d, item = %f\n", i, j, returned_mat[i][j]);
        }
        PyList_SetItem(po_return_mat, i, po_return_mat_i);
    }
    printf("10\n");
    /* Free alocated memory */
    // for (initialize_i = 0; initialize_i < n; initialize_i++)
    // {
    //     for (initialize_j = 0; initialize_j < n; initialize_j++)
    //     {
    //         printf("10.25 %f\n", primary[initialize_i][initialize_j]);
    //     }
    // }
    // for (initialize_i = 0; initialize_i < n; initialize_i++)
    // {
    //     printf("i: %d\n", initialize_i);
    //     free(primary[initialize_i]);
    // }
    // printf("10.5\n");
    // free(primary);
    printf("11\n");
    if (goal == e_spk)
    {
        for (initialize_i = 0; initialize_i < k; initialize_i++)
        {
            free(mu_arr[initialize_i]);
        }
        free(mu_arr);
    }
    printf("12\n");
    return po_return_mat;
}
