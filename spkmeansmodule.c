#define PY_SSIZE_T_CLEAN
#include "spkmeans.h"
#include <Python.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

enum goal_enum {
    spk = 1,
    wam = 2,
    ddg = 3,
    lnorm = 4,
    jacobi = 5,
    kmeans = 6,
};

PyMODINIT_FUNC PyInit_mykmeanssp(void);
static PyObject* fit(PyObject *self, PyObject *args);

/* This part of the code is taken from the presentations of the course */
static PyMethodDef kmeansMethods[] = {
    {"fit", 
    (PyCFunction) fit,
    METH_VARARGS, 
    PyDoc_STR("The main algorithmical function,
    which receives a goal flag and performs the wanted goal")},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef moduleDef = {
    PyModuleDef_HEAD_INIT,
    "myspkmeans",
    NULL, 
    -1, 
    kmeansMethods
};

PyMODINIT_FUNC PyInit_myspkmeans(void)
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
    PyObject* po_primary_i;
    PyObject* po_mu_arr_i;
    PyObject* po_primary_i_j;
    PyObject* po_mu_arr_i_j;
    PyObject* po_mat_i;

    int n, d, k, max, goal_num;
    double EPSILON;
    PyObject* po_primary;
    double** primary;
    PyObject* po_mu_arr;
    double** mu_arr;
    PyObject* po_mat;

    /* Receive the useful information from the user */
    if (!PyArg_ParseTuple(args, "OiiidiOi", &po_primary, &n, &d, &k, 
    &EPSILON, &max, &po_mu_arr, &goal_num))
    {
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }

    if (goal == 1) goal_enum = spk;
    else if (goal == 2) goal_enum = wam;
    else if (goal == 3) goal_enum = ddg;
    else if (goal == 4) goal_enum = lnorm;
    else if (goal == 5) goal_enum = jacobi;
    else if (goal == 6) goal_enum = kmeans;
    else
    {
        /* Unexpected goal value */
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");}
    }

    /* Transform the given PyObjects to double matrices */
    primary = (double **)malloc(n * sizeof(double *));
    if (!primary)
    {
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }

    for (initialize_i = 0; initialize_i < n; initialize_i++)
    {
        primary[initialize_i] = (double *)malloc(d * sizeof(double));
        if (!primary[initialize_i])
        {
            printf("An Error Has Occurred\n");
            return Py_BuildValue("");
        }
    }

    for (i = 0; i < n; i++)
    {
        po_primary_i = PyList_GetItem(po_primary, i);
        for (j = 0; j < d; j++)
        {
            po_primary_i_j = PyList_GetItem(po_primary_i, j);
            primary[i][j] = (double)PyFloat_AsDouble(po_primary_i_j);
        }
    }

    /* Fix from here, add kmeans */
    if (goal_enum == spk)
    {
        mu_arr = (double **)malloc(k * sizeof(double *));
        if (!mu_arr)
        {
            printf("An Error Has Occurred\n");
            return Py_BuildValue("");
        }

        for (initialize_i = 0; initialize_i < k; initialize_i++)
        {
            mu_arr[initialize_i] = (double *)malloc(d * sizeof(double));
            if (!mu_arr[initialize_i])
            {
                printf("An Error Has Occurred\n");
                return Py_BuildValue("");
            }
        }

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
    returned_mat = execute_goal(primary, n, d, &k, EPSILON, max, mu_arr, goal);

    if (goal_enum == jacobi)
    {
        /* Return  */
    }
    else
    {
        /* Return  */
    }

    /* Transform the centroid matrix to a PyObject (PyList) */
    po_final_centroids = PyList_New(k);
    for (i = 0; i < k; i++)
    {
        po_final_centroids_i = PyList_New(d);
        for (j = 0; j < d; j++)
        {
            PyList_SetItem(po_final_centroids_i, j, 
            Py_BuildValue("d", centroids[i][j]));
        }
        PyList_SetItem(po_final_centroids, i, po_final_centroids_i);
    }

    /* Free alocated memory */
    for (initialize_i = 0; initialize_i < n; initialize_i++)
    {
        free(datapoints[initialize_i]);
    }
    for (initialize_i = 0; initialize_i < k; initialize_i++)
    {
        free(centroids[initialize_i]);
    }
    free(datapoints);
    free(centroids);

    return po_final_centroids;
}
