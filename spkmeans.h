void Kmeans(double **matrix, double **mu, int n, int d, int k, int max_iter, double EPSILON);
double dist(double x1, double *x2, int dim);
/ An inner function used to calculate distances by the Kmeans function
