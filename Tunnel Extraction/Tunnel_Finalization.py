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


# Extract ribosome names 
pro_list = []

with open('/Users/Protein_Code.txt', 'r') as f:
    for lines in f.readlines()[1:]:
        line = lines.strip('\n')  
        data_set = line.split()
        pro_list.append(data_set[0])

# Count the number of tunnels (>1) extracted from each species without specified exit coordinates
tunnel_num = []
for i in range (len(pro_list)):
    folder_path = '/Users/%s/pymol/csv' %str(pro_list[i])
    file_count = count_files_in_folder(folder_path)
    tunnel_num.append(file_count)

# Tunnel refinement

sim = [] # Tunnel similarity score
end = [] # Distance between tunnel exit port and conserved residue of uL23/eL39

correct = [] # Select "correct" tunnel (number) based on similarity 
correct_end = [] # Select "correct" tunnel (number) based on exit-port-distance

correct_rad = [] # Exit port radius of 'correct'
end_rad = [] # Exit port radius of 'correct_end'

# Further evaluation of the tunnel structures is required due to blockages or bad extractions
need_check = []

for i in range (len(pro_list)):
    coordx = []
    coordy = []
    coordz = []
    rad = []
    dup_list = []

    # Extract the coords of conserved residues of uL23/eL39
    with open('/Users/exit_coords_%s.txt' %str(pro_list[i]), 'r+') as f:
        lines = f.readlines()
        x_23 = float(lines[0].strip().split()[0])
        y_23 = float(lines[1].strip().split()[0])
        z_23 = float(lines[2].strip().split()[0])

    # Extract the centerline coords of the tunnel extracted with specified exit coordinates
    with open('/Users/%s_single/csv/tunnel_1.csv' %str(pro_list[i]), 'r+') as f:
        for lines in f.readlines()[1:]:
            line = lines.strip('\n')  
            data_set = line.split(',')
            coordx.append(float("{:.2f}".format(float(data_set[4]))))
            coordy.append(float("{:.2f}".format(float(data_set[5]))))
            coordz.append(float("{:.2f}".format(float(data_set[6]))))
            rad.append(float("{:.2f}".format(float(data_set[3]))))
            
    # Remove overlapping coords           
    for k in range (1,len(coordx)):
        if coordx[k] == coordx[k-1] and coordy[k] == coordy[k-1] and coordz[k] == coordz[k-1]:
            dup_list.append(k)
            
    dup_list.sort(reverse=True)
    for l in dup_list:
        del coordx[l]
        del coordy[l]
        del coordz[l]
        del rad[l]

    # Generate a set of new centerline coords at equal distance intervals
    points = np.column_stack((coordx, coordy, coordz))
    distances = np.sqrt(np.sum(np.diff(points, axis=0)**2, axis=1))
    cumulative_distances = np.insert(np.cumsum(distances), 0, 0)
    normalized_distances = cumulative_distances / cumulative_distances[-1]
    t = np.linspace(0, 1, 100) 
    interp = CubicSpline(normalized_distances, points, axis=0)
    new_centerline_points = interp(t)
    coordx_new = new_centerline_points[:,0]
    coordy_new = new_centerline_points[:,1]
    coordz_new = new_centerline_points[:,2]


    loop_num = tunnel_num[i]
    correct_num = [] # Number of candidate tunnel
    sim_min = [] # Similarity score of each candidate tunnel
    rad_min = [] # Exit port radius of each candidate tunnel
    last_dis_min = [] # Distance between exit port and conserved residue of uL23/eL39 of each candidate tunnel
    
    # Extract the tunnel centerline coordinates of each species without specified exit coordinates
    for loop in range (1,loop_num):
        center_x = []
        center_y = []
        center_z = []
        center_rad = []
        center_dup_list = []
        sim_num = 0
        with open('/Users/%s/csv/tunnel_%s.csv' %(str(pro_list[i]), loop), 'r+') as f:
            for lines in f.readlines()[1:]:
                line = lines.strip('\n')  
                data_set = line.split(',')
                center_x.append(float("{:.2f}".format(float(data_set[4]))))
                center_y.append(float("{:.2f}".format(float(data_set[5]))))
                center_z.append(float("{:.2f}".format(float(data_set[6]))))
                center_rad.append(float("{:.2f}".format(float(data_set[3]))))
        
        # Remove overlapping coords for each tunnel
        for k in range (1,len(center_x)):
            if center_x[k] == center_x[k-1] and center_y[k] == center_y[k-1] and center_z[k] == center_z[k-1]:
                center_dup_list.append(k)
            
        center_dup_list.sort(reverse=True)
        for l in center_dup_list:
            del center_x[l]
            del center_y[l]
            del center_z[l]
            del center_rad[l]
            
        # Generate a set of new centerline coords at equal distance intervals for each tunnel
        center_points = np.column_stack((center_x, center_y, center_z))
        center_distances = np.sqrt(np.sum(np.diff(center_points, axis=0)**2, axis=1))
        center_cumulative_distances = np.insert(np.cumsum(center_distances), 0, 0)
        center_normalized_distances = center_cumulative_distances / center_cumulative_distances[-1]
        center_t = np.linspace(0, 1, 100) 
        center_interp = CubicSpline(center_normalized_distances, center_points, axis=0)
        new_center = center_interp(center_t)
        center_x_new = new_center[:,0]
        center_y_new = new_center[:,1]
        center_z_new = new_center[:,2]

        # Compute the similarity score between this tunnel and the tunnel derived from the specified exit coordinates 
        # by calculating the deviations in the coordinates of each centerline
        for j in range (len(center_x_new)):
            devi = np.sqrt((center_x_new[j] - coordx_new[j])**2+(center_y_new[j] - coordy_new[j])**2+(center_z_new[j] - coordz_new[j])**2)
            sim_num += devi
                    
        similarity = sim_num 
        sim_min.append(similarity)
        
        correct_num.append(loop)
        
        rad_min.append(center_rad[-1])
        
        last_x = center_x_new[-1]
        last_y = center_y_new[-1]
        last_z = center_z_new[-1]

        # Compute the Euclidean distance between the exit port and the conserved residue
        last_dis = np.sqrt((last_x - x_23)**2+(last_y - y_23)**2+(last_z - z_23)**2)
        last_dis_min.append(last_dis)
   
    # Find the two candidate
    smallest, indices = find_two_smallest(sim_min)
    
    last_dis_smallest, dis_indices = find_two_smallest(last_dis_min)
   
    sim.append(smallest)
    end.append(last_dis_smallest)
    
    correct.append([correct_num[indices[0]],correct_num[indices[1]],correct_num[indices[2]]])
    correct_end.append([correct_num[dis_indices[0]],correct_num[dis_indices[1]],correct_num[dis_indices[2]]])
    
    correct_rad.append([rad_min[indices[0]],rad_min[indices[1]],rad_min[indices[2]]])
    end_rad.append([rad_min[dis_indices[0]],rad_min[dis_indices[1]],rad_min[dis_indices[2]]])
    
    if correct_num[indices[0]] != correct_num[dis_indices[0]]:
        need_check.append(i)



