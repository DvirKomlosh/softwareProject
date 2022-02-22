import sys
import numpy as np
from numpy import linalg as lng
from enum import Enum
import myspkmeans

np.random.seed(0)

# The requested goal enum
class goal_enum(Enum):
    spk = 1
    wam = 2
    ddg = 3
    lnorm = 4
    jacobi = 5


# This function manages the reading from the input file,
# and creation of the data matrix.
def input_file_management(file_name):
    # Open the given file
    input_file = None
    try:
        input_file = open(file_name, "r")
    except FileNotFoundError:
        print("An Error Has Occurred")
        exit(-1)

    lines = input_file.readlines()
    input_file.close()
    N = len(lines)  # The amount of data points
    # Set dataset to have the size of the final data matrix
    dataset = np.zeros((N, len(lines[0].split(','))))  
    
    # Fill up the dataset with the file's data
    for row_index in range(N):
        curr_line_vector = lines[row_index].strip('\n').split(',')
        for column_index in range(len(curr_line_vector)):
            dataset[row_index][column_index] = curr_line_vector[column_index]

    return dataset


# This function preforms the kmeans++ algorithm.
# It receives the amount of clusters and the data,
# and returns the centroids indexes, and the centroids themselves.
def k_means_pp(k, clustering_data):
    n = len(clustering_data)
    d = np.zeros(n)
    p = np.ones(n)
    mu = np.zeros(shape=(k, len(clustering_data[0])))
    mu_indecies = np.zeros(k)

    mu_indecies[0] = np.random.choice(a=np.arange(n))
    mu[0] = clustering_data[int(mu_indecies[0])]

    for index in range(1, k):

        for l_index in range(n):
            d[l_index] = np.min([lng.norm(clustering_data[l_index] - \
                mu[j_index]) ** 2 for j_index in range(index)])

        d_sum = sum(d)
        for l_index in range(n):
            p[l_index] = d[l_index] / d_sum

        mu_indecies[index] = np.random.choice(a=np.arange(n), p=p)
        mu[index] = clustering_data[int(mu_indecies[index])]

    return mu_indecies, mu


# This function is responsible of outputting the data of the spk procedure.
def output_spk(indices, centroids):
    print(",".join(["{:.4f}".format(indices[i]) for i in range(len(indices))]))
    for row_index in range(len(centroids)):
        print(",".join(["{:.4f}".format(centroids[row_index][i]) \
            for i in range(len(centroids[row_index]))]))


# This function is responsible of outputting the data of all other procedures.
def output_other_than_spk(mat):
    for row_index in range(len(mat)):
        print(",".join(["{:.4f}".format(mat[row_index][i]) \
            for i in range(len(mat[row_index]))]))


if __name__ == '__main__':
    # Obtain parameters from command
    params_from_cmd = list(sys.argv)
    if len(params_from_cmd) != 4:
        # Invalid amount of parameters
        print("Invalid Input!")
        exit(1)

    max_rotations = 100
    max_iterations = 300
    jacobi_epsilon = 10 ** -15
    kmeans_epsilon = 0

    # Processing of the first parameter - k
    k = -1
    if params_from_cmd[1].isnumeric():
        k = int(params_from_cmd[1])
    else:
        # Invalid type of the first parameter
        print("Invalid Input!")
        exit(1)
    if k < 0:
        # Invalid value of the first parameter
        print("Invalid Input!")
        exit(1)

    # Processing of the second parameter - goal
    goal_string = params_from_cmd[2]
    if goal_string == "spk":
        goal = goal_enum.spk
    elif goal_string == "wam":
        goal = goal_enum.wam
    elif goal_string == "ddg":
        goal = goal_enum.ddg
    elif goal_string == "lnorm":
        goal = goal_enum.lnorm
    elif goal_string == "jacobi":
        goal = goal_enum.jacobi
    else:
        # Invalid value of the second parameter
        print("Invalid Input!")
        exit(1)

    # Processing of the third parameter - file_name
    file_name = params_from_cmd[3]
    #  The file extension is .txt or .csv
    if file_name.split(".")[-1] not in ["txt", "csv"]:
        # Invalid value of the third parameter
        print("Invalid Input!")
        exit(1)
    data = input_file_management(file_name)
    N = len(data)
    d = len(data[0])

    if k >= N:
        # Invalid value of the first parameter - the number of clusters 
        # should be smaller then the number of data points.
        print("Invalid Input!")
        exit(1)

    # Activate the wanted goal function using the CAPI file.
    if goal == goal_enum.spk:
        # Perform full spectral kmeans
        jacobi_mat = myspkmeans.jacobi(data, k, jacobi_epsilon, max_rotations)
        mu_indices, mu = k_means_pp(k, jacobi_mat)
        found_centroids = myspkmeans.fit(N, d, k, max_iterations,
            kmeans_epsilon, np.ndarray.tolist(data), np.ndarray.tolist(mu))
        if not found_centroids:
            exit(1)
        output_spk(mu_indices, found_centroids)
    else:
        if goal == goal_enum.wam:
            # Calculate and output the Weighted Adjacency Matrix
            mat = myspkmeans.wam(data)
        elif goal == goal_enum.ddg:
            # Calculate and output the Diagonal Degree Matrix
            mat = myspkmeans.ddg(data)
        elif goal == goal_enum.lnorm:
            # Calculate and output the Normalized Graph Laplacian
            mat = myspkmeans.lnorm(data)
        elif goal == goal_enum.jacobi:
            # Calculate and output the eigenvalues and eigenvectors
            mat = myspkmeans.jacobi(data, k, jacobi_epsilon, max_rotations)
        else:
            # Invalid goal value
            print("An Error Has Occurred")
            exit(1)
        output_other_than_spk(mat)
