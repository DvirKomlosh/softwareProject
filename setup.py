import sys
import numpy as np
from numpy import linalg as lng
# Hello

def input_file_management(args, first_index):
    # Open the given file
    input_f_1 = None
    try:
        input_f_1 = open(args[first_index], "r")
    except FileNotFoundError:
        print("An Error Has Occurred")
        exit(-1)

    lines = input_f_1.readlines()
    length = len(lines)  # check the amounts of data points
    dataset = np.zeros(
        (length, 2 * len(lines[0].split(',')) - 2))  # set data set to have the size of the final data matrix

    input_f_1.close()
    input_f_1 = open(args[first_index], "r")

    input_f_2 = None
    try:
        input_f_2 = open(args[first_index + 1], "r")
    except FileNotFoundError:
        print("An Error Has Occurred")
        exit(-1)

    # Fill up the dataset with the file's data
    finished_reading = False
    while not finished_reading:
        # Read a line from the given file
        curr_line = input_f_1.readline()
        if curr_line == '':
            # The program read the last line of the file
            finished_reading = True
            continue
        # The program read the last line of the file
        line_vector = curr_line.strip('\n').split(',')
        for index_j in range(len(line_vector) - 1):
            dataset[int(float(line_vector[0]))][index_j] = float(
                line_vector[index_j + 1])  # placing the half-vector into the right positions in the matrix
    input_f_1.close()

    finished_reading = False
    while not finished_reading:
        # Read a line from the given file
        curr_line = input_f_2.readline()
        if curr_line == '':
            # The program read the last line of the file
            finished_reading = True
            continue
        line_vector = curr_line.strip('\n').split(',')
        for index_j in range(len(line_vector) - 1):
            dataset[int(float(line_vector[0]))][int(index_j + (len(line_vector) - 1))] = float(
                line_vector[index_j + 1])  # placing the half-vector into the right positions in the matrix

    input_f_2.close()

    return np.array(dataset)


def k_means_pp(k, clustering_data):
    n = len(clustering_data)
    d = np.zeros(n)
    p = np.ones(n)
    mu = np.zeros(shape=(k, len(clustering_data[0])))
    mu_indecies = np.zeros(k)

    np.random.seed(0)

    mu_indecies[0] = np.random.choice(a=np.arange(n))
    mu[0] = clustering_data[int(mu_indecies[0])]

    for index in range(1, k):

        for l_index in range(n):
            d[l_index] = np.min([lng.norm(clustering_data[l_index] - mu[j_index]) ** 2 for j_index in range(index)])

        d_sum = sum(d)
        for l_index in range(n):
            p[l_index] = d[l_index] / d_sum

        mu_indecies[index] = np.random.choice(a=np.arange(n), p=p)
        mu[index] = clustering_data[int(mu_indecies[index])]

    return mu_indecies, mu


if __name__ == '__main__':
    # Obtain parameters from command
    params_from_cmd = list(sys.argv)
    if len(params_from_cmd) != 5 and len(params_from_cmd) != 6:
        # Invalid amount of parameters
        print("Invalid Input!")
        exit(1)

    # Processing of the first two parameters
    amount_of_clusters = -1
    max_iteration = -1
    input_file_param_index = 2

    if params_from_cmd[1].isnumeric():
        amount_of_clusters = int(params_from_cmd[1])
    else:
        # Invalid type of the first parameter
        print("Invalid Input!")
        exit(1)

    max_iteration = 300
    if len(params_from_cmd) == 6:
        if params_from_cmd[2].isnumeric():
            max_iteration = int(params_from_cmd[2])
            input_file_param_index = 3
        else:
            # Invalid type of the first parameter
            print("Invalid Input!")
            exit(1)
    epsilon = 0
    try:
        epsilon = float(params_from_cmd[input_file_param_index])
    except ValueError:
        # Invalid type of the first parameter
        print("Invalid Input!")
        exit(1)

    if amount_of_clusters < 1 or max_iteration < 1 or epsilon < 0:
        # Invalid value of parameter
        print("Invalid Input!")
        exit(1)

    # Processing of the third/fourth parameter
    data = input_file_management(params_from_cmd, input_file_param_index + 1)

    if len(data) < amount_of_clusters:
        # Invalid value of the first parameter, the number of clusters should be smaller then the number of data points.
        print("An Error Has Occurred")
        exit(1)

    # Reformat the centroids data to the wanted format
    initialized_centroids_indecies, initialized_centroids = k_means_pp(amount_of_clusters, data)

    found_centroids = mykmeanssp.fit(len(data), len(data[0]), amount_of_clusters, max_iteration,
                                     epsilon, np.ndarray.tolist(data), np.ndarray.tolist(initialized_centroids))
    if not found_centroids:
        exit(1)

    print(','.join([str(int(x)) for x in initialized_centroids_indecies]))

    for i in range(len(found_centroids)):
        for j in range(len(found_centroids[i])):
            found_centroids[i][j] = "{:.4f}".format(found_centroids[i][j])

    for i in range(len(found_centroids)):
        print(','.join(found_centroids[i]))
