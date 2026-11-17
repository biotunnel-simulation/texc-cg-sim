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

def find_three_smallest(numbers):
    sorted_numbers = sorted(numbers)
    smallest_three = sorted_numbers[:3]
    index_of_smallest = [numbers.index(num) for num in smallest_three]
    return smallest_three, index_of_smallest


# Extract ribosome PDB IDs 
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


# Tunnel refinement process
sim = [] # Tunnel dissimilarity score
end = [] # Distance between tunnel exit port and conserved residue of uL23/eL39

correct = [] # "correct" tunnel (number) based on dissimilarity 
correct_end = [] # "correct" tunnel (number) based on exit-port distance

correct_rad = [] # Exit-port radius of 'correct'
end_rad = [] # Exit-port radius of 'correct_end'

# Downstream analysis is required for the following tunnels
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

        # Compute the dissimilarity score between this tunnel and the tunnel derived from the specified exit coordinates 
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
   
    # Identify the three smallest dissimilarity scores and the tunnel number
    smallest, indices = find_three_smallest(sim_min)

    # Identify the three smallest exit-port distances and the tunnel number
    last_dis_smallest, dis_indices = find_three_smallest(last_dis_min)
   
    sim.append(smallest)
    end.append(last_dis_smallest)

    # Tunnel numbers based on smallest dissimilarity scores
    correct.append([correct_num[indices[0]],correct_num[indices[1]],correct_num[indices[2]]])

    # Tunnel numbers based on smallest exit-port distances
    correct_end.append([correct_num[dis_indices[0]],correct_num[dis_indices[1]],correct_num[dis_indices[2]]])

    # Exit-port radius of 'correct'
    correct_rad.append([rad_min[indices[0]],rad_min[indices[1]],rad_min[indices[2]]])

    # Exit-port radius of 'correct_end'
    end_rad.append([rad_min[dis_indices[0]],rad_min[dis_indices[1]],rad_min[dis_indices[2]]])

    # The “ideal” candidate tunnel is assumed to have both the smallest dissimilarity scores and exit-port distances
    # If no “ideal” candidate tunnel exists, perform a downstream analysis on these tunnels
    if correct_num[indices[0]] != correct_num[dis_indices[0]]:
        need_check.append(i)


# Downstream analysis
# The sum of n, nn, and nnn should equal the total number of species
nn = 0
nnn = 0
nnnn = 0
final_tunnel = []

for i in range (len(pro_list)):
    # If the tunnel number with the smallest dissimilarity scores = the tunnel number with the smallest exit-port distance
    if correct[i][0] == correct_end[i][0]:
        
        # Then check whether the “ideal” candidate tunnel has the “ideal” exit-port radius, 6 was found to be a good standard
        if correct_rad[i][0] < 6: 
            
            # If the exit port radius < 6 angstrom
            # The final tunnel assumed to be the one with the samllest dissimilarity scores for now
            output_value = correct[i][0]
            
            # Then further check whether the tunnel number associated with the second-lowest dissimilarity score appears in the list of tunnel numbers ranked by the smallest exit-port distances
            if correct[i][1] in correct_end[i]:
                
                # Then check whether the tunnel number associated with the third-lowest dissimilarity score appears in the list of tunnel numbers ranked by the smallest exit-port distances
                if correct[i][2] in correct_end[i]:
                    
                    # If the tunnel number associated with the third-lowest dissimilarity score also appears in the list of tunnel numbers ranked by the smallest exit-por
                    # The final tunnel considered the one with the largest exit-port radius
                    
                    maximum = max(correct_rad[i])
                    max_index = correct_rad[i].index(maximum) 
                    output_value = correct[i][max_index]
                    
                else:

                    # If not, only consider the tunnels with the samllest and the second-lowest dissimilarity score, then further check their exit-port radius
                    if correct_rad[i][1] >= correct_rad[i][0]:

                        # In this case, the final tunnel considered the one with the larger exit-port radius
                        output_value = correct[i][1]

                # Otherwise, the tunnel with the smallest dissimilarity score is selected as the final tunnel
                final_tunnel.append(int(output_value)) 
                    
            # If only the tunnel number associated with the third-lowest dissimilarity score appears in the list of tunnel numbers ranked by the smallest exit-port distances
            elif correct[i][2] in correct_end[i]:
                
                # Then compare their exit-port radius between tunnels with the samllest and the third-lowest dissimilarity score
                # The final tunnel considered the one with the larger exit-port radius
                if correct_rad[i][2] >= correct_rad[i][0]:
                    output_value = correct[i][2]
                    
                final_tunnel.append(int(output_value)) 
                
            else:
                final_tunnel.append(int(output_value)) 
     
        else:
            # The tunnel is considered the final tunnel if its exit-port radius is >= 6 angstrom
            final_tunnel.append(int(correct[i][0]))
        nn += 1  

    # If the tunnel number with the smallest dissimilarity score differs from that with the smallest exit-port distance
    # Check whether the first two tunnel numbers are swapped between the two lists
    elif correct[i][0] == correct_end[i][1] and correct[i][1] == correct_end[i][0]:

        # If so, select the tunnel with the larger exit-port radius as the final tunnel
        if correct_rad[i][0] >= correct_rad[i][1]:
            final_tunnel.append(int(correct[i][0]))
        if correct_rad[i][0] < correct_rad[i][1]:
            final_tunnel.append(int(correct[i][1]))

        nnn += 1 

    # If the first two tunnel numbers are not swapped (not the same numbers)
    else:

        # Identify the tunnel number common to both lists
        common_numbers = set(correct[i]) & set(correct_end[i]) 
        common_indices = [j for j, num in enumerate(correct[i]) if num in common_numbers]

        # If only one tunnel number is common to both lists, select it as the final tunnel
        if len(common_indices) == 1:
            final_tunnel.append(int(correct[i][common_indices[0]]))
            
        # If two tunnel numbers appears in both lists
        if len(common_indices) == 2:

            # Identify the one with the larger exit-port radius
            maximum = max(correct_rad[i][common_indices[0]], correct_rad[i][common_indices[1]])

            # If its exit-port radius <= 6 angstrom, select the tunnel with the largest exit-port radius from the two lists
            if maximum <= 6:
                combined_list = correct[i] + correct_end[i]
                combined_list2 = correct_rad[i] + end_rad[i]
                maximum = max(combined_list2)
                max_index = combined_list2.index(maximum) 
                output_value = combined_list[max_index]
                final_tunnel.append(int(output_value)) 

            # If its exit-port radius > 6 angstrom, it is considered as the final tunnel
            else:  
                max_index = correct_rad[i].index(maximum) 
                output_value = correct[i][max_index]
                final_tunnel.append(int(output_value))
            
        # If all tunnel numbers appears in both lists but in different order
        # The final tunnel considered the one with the largest exit-port radius
        if len(common_indices) == 3:
            maximum = max(correct_rad[i][common_indices[0]], correct_rad[i][common_indices[1]], correct_rad[i][common_indices[2]])
            max_index = correct_rad[i].index(maximum) 
            output_value = correct[i][max_index]
            final_tunnel.append(int(output_value))

        # A very rare and exceptional case: there is no common tunnel number between the two lists
        # Select the tunnel with the largest exit-port radius as the final tunnel
        # Note: Visual inspection may be required !!!!!
        if len(common_indices) == 0:
            combined_list = correct[i] + correct_end[i]
            combined_list2 = correct_rad[i] + end_rad[i]
            maximum = max(combined_list2)
            max_index = combined_list2.index(maximum) 
            output_value = combined_list[max_index]
            final_tunnel.append(int(output_value)) 
        nnnn += 1 


# Once the final tunnel is selected for each species, the 2D radial plot of each tunnel can be generated
# A 2D plot can serve as an additional validation method to assess the accuracy of the previous tunnel selection process and improve the robustness of the results
for i in range (len(pro_list)):
    center_x = []
    center_y = []
    center_z = []
    center_rad = []
    center_dup_list = []
    
    with open('/Users/%s/csv/tunnel_%s.csv' %(str(pro_list[i]), final_tunnel[i]), 'r+') as f:
        for lines in f.readlines()[1:]:
            line = lines.strip('\n')  
            data_set = line.split(',')
            center_x.append(float("{:.2f}".format(float(data_set[4]))))
            center_y.append(float("{:.2f}".format(float(data_set[5]))))
            center_z.append(float("{:.2f}".format(float(data_set[6]))))
            center_rad.append(float("{:.2f}".format(float(data_set[3]))))
            if float("{:.2f}".format(float(data_set[3])))>10.5:
                print(str(pro_list[i]))
                
    # Remove overlapping coords         
    for k in range (1,len(center_x)):
        if center_x[k] == center_x[k-1] and center_y[k] == center_y[k-1] and center_z[k] == center_z[k-1]:
            center_dup_list.append(k)     
    center_dup_list.sort(reverse=True)
    for l in center_dup_list:
        del center_x[l]
        del center_y[l]
        del center_z[l]
        del center_rad[l]
        
    coord = [0] 
    sum_cor = 0
    
    for v in range (1,len(center_x)):
        delt = np.sqrt((center_x[v]-center_x[v-1])**2+(center_y[v]-center_y[v-1])**2+(center_z[v]-center_z[v-1])**2)
        sum_cor += delt 
        coord.append(sum_cor)

    # Generate the 2D plot
    plt.plot(coord, center_rad,linewidth = 1)
    plt.xlabel('Tunnel length')
    plt.ylabel('Radius')
