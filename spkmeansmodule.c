#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
    PyObject* po_additional_i;
    PyObject* po_primary_i_j;
    PyObject* po_additional_i_j;
    PyObject* po_mat_i;

    int n, d, k, max, goal;
    double EPSILON;
    PyObject* po_primary;
    double** primary;
    PyObject* po_additional;
    double** additional;
    PyObject* po_mat;

    /* Receive the useful information from the user */
    if (!PyArg_ParseTuple(args, "OiiidiOi", &po_primary, &n, &d, &k, 
    &EPSILON, &max, &po_additional, &goal))
    {
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }

    /* Transform the given PyObjects to double matrices */
    primary = (double **)malloc(n * sizeof(double *)); 
    additional = (double **)malloc(k * sizeof(double *));
    if (!primary || !additional)
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
    for (initialize_i = 0; initialize_i < k; initialize_i++)
    {
        additional[initialize_i] = (double *)malloc(d * sizeof(double));
        if (!additional[initialize_i])
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

    for (i = 0; i < k; i++)
    {
        po_additional_i = PyList_GetItem(po_additional, i);
        for (j = 0; j < d; j++)
        {
            po_additional_i_j = PyList_GetItem(po_additional_i, j);
            additional[i][j] = (double)PyFloat_AsDouble(po_additional_i_j);
        }
    }

    /* Activate the correct goal function */

    /* Activate the KMeans() function */
    Kmeans(datapoints, centroids, n, d, k, max_iter, EPSILON);

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




