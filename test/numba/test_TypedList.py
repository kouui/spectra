import numpy as np
from numba import njit
from numba.typed import List
import time

def make_array_lists(shapes):

    typed_list = List()
    reflected_list = []

    for shape in shapes:
        typed_list.append(np.zeros(shape, dtype=np.float64))
        reflected_list.append(np.zeros(shape, dtype=np.float64))

    return typed_list, reflected_list

@njit(cache=True)
def mutate_array_list(array_list):

    for j in range(len(array_list)):
        dim0, dim1 = array_list[j].shape
        for row in range(dim0):
            for col in range(dim1):
                array_list[j][row, col] += 1


if __name__ == "__main__":

    array_shapes = [(200,200), (300,300), (400,400)]
    nrun = 1000

    typed_list, reflected_list = make_array_lists(array_shapes)

    t0 = time.time()
    [mutate_array_list(typed_list) for i in range(nrun)]
    average_time = (time.time() - t0)/nrun
    print("Typed list took {:.6f} seconds per run.".format(average_time))

    t0 = time.time()
    [mutate_array_list(reflected_list) for i in range(nrun)]
    average_time = (time.time() - t0)/nrun
    print("Reflected list took {:.6f} seconds per run.".format(average_time))
