import os
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import numpy as np
import statistics
import shutil

def count_files_in_folder(folder_path):
    count = 0
    for _, _, files in os.walk(folder_path):
        count += len(files)
    return count

def find_two_smallest(numbers):
    sorted_numbers = sorted(numbers)
    smallest_two = sorted_numbers[:3]
    index_of_smallest = [numbers.index(num) for num in smallest_two]
    return smallest_two, index_of_smallest

def make_lists_equal_length(list1, list2, list3, list4):
    len_diff = abs(len(list1) - len(list2))
    if len(list1) < len(list2):
        list1 += [0] * len_diff
    else:
        list2 += [0] * len_diff
        list3 += [0] * len_diff
        list4 += [0] * len_diff
    return list1, list2, list3, list4


pro_list = []

with open('/Users/Protein_Code.txt', 'r') as f:
    for lines in f.readlines()[1:]:
        line = lines.strip('\n')  
        data_set = line.split()
        pro_list.append(data_set[0])


