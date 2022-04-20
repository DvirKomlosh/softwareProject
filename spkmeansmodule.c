#include "spkmeans.h"
#include "algorithm.c"

#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject* fit(PyObject *self, PyObject *args);

/* This part of the code is taken mostly from the course's presentations */
static PyMethodDef kmeansMethods[] = {
    {"fit", 
    (PyCFunction) fit,
    METH_VARARGS, 
    PyDoc_STR("The main algorithmical function, 
    which receives a goal flag and other required parameters 
    and performs the wanted goal")},
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
    
    /* Indexers */
    int initialize_i;
    int i, j;

    /* Variables which help the allocation/filling process */
    PyObject* po_primary_i;
    PyObject* po_mu_arr_i;
    PyObject* po_primary_i_j;
    PyObject* po_mu_arr_i_j;
    PyObject* po_return_mat_i;

    enum goal_enum goal; /* The chosen goal */
    int n, d, k, goal_num; /* Various parameters from the Python program */

    PyObject* po_primary; /* The PyObject of the primary dataset */
    double** primary; /* The primary dataset */
    PyObject* po_mu_arr; /* The PyObject of the centroids */
    double** mu_arr = NULL; /* The centroids */
    
    /* The PyObject of the matrix, returned from the execute_goal() function */
    PyObject* po_return_mat; 
    /* The matrix returned from the execute_goal() function */
    double** returned_mat;  

    /* Receive the information from the user */
    if (!PyArg_ParseTuple(args, "OiiiOi", &po_primary, &n, &d, &k, 
    &po_mu_arr, &goal_num))
    {
        /* An error has occurred during the information extraction */
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }

    if (goal_num == 1) goal = e_spk;
    else if (goal_num == 2) goal = e_wam;
    else if (goal_num == 3) goal = e_ddg;
    else if (goal_num == 4) goal = e_lnorm;
    else if (goal_num == 5) goal = e_jacobi;
    else if (goal_num == 6) goal = e_kmeans;
    else
    {
        /* An unexpected goal value is detected */
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }

    /* Parse the primary matrix */
    
    /* Allocate the matrix */
    primary = (double **)malloc(n * sizeof(double *));
    if (!primary)
    {
        /* The allocation has failed */
        printf("An Error Has Occurred\n");
        return Py_BuildValue("");
    }

    for (initialize_i = 0; initialize_i < n; initialize_i++)
    {
        /* Allocate each inner array */
        primary[initialize_i] = (double *)malloc(d * sizeof(double));
        if (!primary[initialize_i])
        {
            /* The allocation has failed */
            printf("An Error Has Occurred\n");
            return Py_BuildValue("");
        }
    }
    
    /* Fill the data from the corresponding PyObject */ 
    for (i = 0; i < n; i++)
    {
        po_primary_i = PyList_GetItem(po_primary, i);
        for (j = 0; j < d; j++)
        {
            po_primary_i_j = PyList_GetItem(po_primary_i, j);
            primary[i][j] = (double)PyFloat_AsDouble(po_primary_i_j);
        }
    }

    if (goal == e_spk)
    {   
        /* Parse the centroids array (which is a matrix) */

        /* Allocate the matrix */
        mu_arr = (double **)malloc(k * sizeof(double *));
        if (!mu_arr)
        {
            /* The allocation has failed */
            printf("An Error Has Occurred\n");
            return Py_BuildValue("");
        }

        for (initialize_i = 0; initialize_i < k; initialize_i++)
        {
            mu_arr[initialize_i] = (double *)malloc(d * sizeof(double));
            if (!mu_arr[initialize_i])
            {
                /* The allocation has failed */
                printf("An Error Has Occurred\n");
                return Py_BuildValue("");
            }
        }

        /* Fill the data from the corresponding PyObject */ 
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

    /* Transform the returned matrix to a PyObject (PyList) */   
    if (goal == e_spk) po_return_mat = PyList_New(k);
    else if (goal == e_jacobi) po_return_mat = PyList_New(n+1);
    else po_return_mat = PyList_New(n);
    for (i = 0; i < PyList_Size(po_return_mat); i++)
    {
        if (goal == e_spk) po_return_mat_i = PyList_New(d);
        else if (goal == e_kmeans) po_return_mat_i = PyList_New(k);
        else po_return_mat_i = PyList_New(n);
        for (j = 0; j < PyList_Size(po_return_mat_i); j++)
        {
            /* Fill the data */
            PyList_SetItem(po_return_mat_i, j,
                Py_BuildValue("d", returned_mat[i][j]));
        }
        PyList_SetItem(po_return_mat, i, po_return_mat_i);
    }

    /* Free the allocated memory */
    free_matrix(primary, n);
    if (goal == e_spk) free_matrix(mu_arr, k);

    return po_return_mat;
}
