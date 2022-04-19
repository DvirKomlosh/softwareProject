import sys
import numpy as np
from numpy import linalg as lng
from enum import Enum
import spkmeans

np.random.seed(0)  # Set np.random.seed(0) at the beginning of your code.


# The requested goal enum
class GoalEnum(Enum):
    spk = 1
    wam = 2
    ddg = 3
    lnorm = 4
    jacobi = 5
    kmeans = 6


# This function manages the reading from the input file,
# and creation of the data matrix.
def input_file_management(filename):
    # Open the given file
    input_file = None
    try:
        input_file = open(filename, "r")
    except FileNotFoundError:
        print("An Error Has Occurred")
        exit(1)

    lines = input_file.readlines()
    input_file.close()
    N_val = len(lines)  # The amount of data points
    # Set dataset to have the size of the final data matrix
    dataset = np.zeros((N_val, len(lines[0].split(','))))

    # Fill up the dataset with the file's data
    for row_index in range(N_val):
        # Try to transform the input strings into float numbers.
        curr_line_vector = None
        try:
            curr_line_vector = [float(num) for num in lines[row_index].strip('\n').split(',')]
        except ValueError:
            print("Invalid Input! (1)")
            exit(1)
        for column_index in range(len(curr_line_vector)):
            dataset[row_index][column_index] = curr_line_vector[column_index]
    return dataset


# This function preforms the kmeans++ algorithm.
# It receives the amount of clusters and the data,
# and returns the centroids indexes, and the centroids themselves.
def k_means_pp(k_val, clustering_data):
    n_val = len(clustering_data)
    d_val = np.zeros(n_val)
    p = np.ones(n_val)
    mu_arr = np.zeros(shape=(k_val, len(clustering_data[0])))
    indices = np.zeros(k_val)

    indices[0] = np.random.choice(a=np.arange(n_val))
    mu_arr[0] = clustering_data[int(indices[0])]

    for index in range(1, k_val):

        for l_index in range(n_val):
            d_val[l_index] = \
                np.min([lng.norm(clustering_data[l_index] - mu_arr[j_index])
                        ** 2 for j_index in range(index)])

        d_sum = sum(d_val)
        for l_index in range(n_val):
            p[l_index] = d_val[l_index] / d_sum

        indices[index] = np.random.choice(a=np.arange(n_val), p=p)
        mu_arr[index] = clustering_data[int(indices[index])]

    return indices, mu_arr


# This function is responsible of outputting the data of the spk procedure.
def output_spk(indices, centroids):
    print(",".join(["{:.4f}".format(indices[i]) for i in range(len(indices))]))
    for row_index in range(len(centroids)):
        print(",".join(["{:.4f}".format(centroids[row_index][i])
                        for i in range(len(centroids[row_index]))]))


# This function is responsible of outputting the data of all other procedures.
def output_other_than_spk(matrix):
    for row_index in range(len(matrix)):
        print(",".join(["{:.4f}".format(matrix[row_index][i])
                        for i in range(len(matrix[row_index]))]))


# This function receives a matrix, and checks if it is a symmetric matrix.
def jacobi_input_validation(input_mat):
    for i in range(len(input_mat)):
        for j in range(i + 1, len(input_mat)):
            # An asymmetry is detected.
            if input_mat[i][j] != input_mat[j][i]:
                return False
    return True


if __name__ == '__main__':
    # Obtain parameters from command
    params_from_cmd = list(sys.argv)
    if len(params_from_cmd) != 4:
        # Invalid amount of parameters
        print("Invalid Input! (2)")
        exit(1)

    # Pass manually and not as parameters!
    max_rotations = 100
    max_iterations = 300
    jacobi_epsilon = 10 ** -5
    kmeans_epsilon = 0

    # Processing of the second parameter - goal
    goal_string = params_from_cmd[2]
    goal = None
    if goal_string == "spk":
        goal = GoalEnum.spk
    elif goal_string == "wam":
        goal = GoalEnum.wam
    elif goal_string == "ddg":
        goal = GoalEnum.ddg
    elif goal_string == "lnorm":
        goal = GoalEnum.lnorm
    elif goal_string == "jacobi":
        goal = GoalEnum.jacobi
    else:
        # Invalid value of the second parameter
        print("Invalid Input! (5)")
        exit(1)

    # Processing of the third parameter - file_name
    file_name = params_from_cmd[3]
    #  The file extension is .txt or .csv
    if file_name.split(".")[-1] not in ["txt", "csv"]:
        # Invalid value of the third parameter
        print("Invalid Input! (6)")
        exit(1)
    data = input_file_management(file_name)
    N = len(data)
    d = len(data[0])

    if N == 0 or d == 0:
        # The input file is empty
        print("Invalid Input! (7)")
        exit(1)

    # Activate the wanted goal function using the CAPI file.
    if goal == GoalEnum.spk:
        # Processing of the first parameter - k
        k = -1
        if params_from_cmd[1].isnumeric():
            k = int(params_from_cmd[1])
        else:
            # Invalid type of the first parameter
            print("Invalid Input! (3)")
            exit(1)
        if k < 0:
            # Invalid value of the first parameter
            print("Invalid Input! (4)")
            exit(1)
        if k >= N and goal == GoalEnum.spk:
            # Invalid value of the first parameter - the number of clusters 
            # should be smaller then the number of data points.
            print("Invalid Input! (8)")
            exit(1)
        # Perform full spectral kmeans
        T = spkmeans.fit(np.ndarray.tolist(data), 
            N, d, k, None, GoalEnum.kmeans.value)
        mu_indices, mu = k_means_pp(k, T)
        found_centroids = spkmeans.fit(np.ndarray.tolist(data), 
            N, d, k, np.ndarray.tolist(mu), GoalEnum.spk.value)
        if not found_centroids:
            exit(1)
        output_spk(mu_indices, found_centroids)
    else:
        mat = None
        if goal == GoalEnum.wam:
            # Calculate and output the Weighted Adjacency Matrix
            mat = spkmeans.fit(np.ndarray.tolist(data), 
                N, d, -1, None, GoalEnum.wam.value)
        elif goal == GoalEnum.ddg:
            # Calculate and output the Diagonal Degree Matrix
            mat = spkmeans.fit(np.ndarray.tolist(data),
                N, d, -1, None, GoalEnum.ddg.value)
        elif goal == GoalEnum.lnorm:
            # Calculate and output the Normalized Graph Laplacian
            mat = spkmeans.fit(np.ndarray.tolist(data),
                N, d, -1, None, GoalEnum.lnorm.value)
        elif goal == GoalEnum.jacobi:
            # Calculate and output the eigenvalues and eigenvectors
            if not jacobi_input_validation(data):
                # Invalid jacobi matrix
                print("Invalid Input! (9)")
                exit(1)
            mat = spkmeans.fit(np.ndarray.tolist(data), 
                N, d, -1, None, GoalEnum.jacobi.value)
        else:
            # Invalid goal value
            print("An Error Has Occurred")
            exit(1)
        output_other_than_spk(mat)
