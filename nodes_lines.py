import math
import numpy as np 

thetas = np.array([i for i in range(0, 315 + 1 + 1, 45)])
thetas = np.deg2rad(thetas)
names = [f"Angle {theta} degrees" for theta in np.degrees(thetas)]


h = 15
H = 20
D = 40 

def Oct(R,a,z):
    #R= radius of circle, a=angle (0 or 22.5), z = heigth
    layer_nodes = []
    for i in range(8):
        angle = math.radians(a + 45 * i)
        x = R * math.cos(angle)
        y = R * math.sin(angle)
        layer_nodes.append((x, y, z))
    return layer_nodes


all_nodes = [] 
#Define the nodes in each heigth
all_nodes.append(Oct(D/2,22.5,0)) #1
all_nodes.append(Oct(D/2,22.5,h)) #2
all_nodes.append(Oct(D/2,22.5,h+H)) #3
all_nodes.append(Oct(D/2,22.5, 2*h + H)) #4
all_nodes.append(Oct(D/2,22.5, 2*h + 2*H)) #5
all_nodes.append(Oct(D/2,22.5, 3*h + 2*H)) #6




all_lines = []

def create_lines2(nodes, layer_indices):

    num_layers = len(layer_indices)

    for i in range(num_layers - 1):  # Iterate through specified layers
        current_layer = nodes[layer_indices[i]]
        next_layer = nodes[layer_indices[i + 1]]

        for j in range(len(current_layer)):
            current_point = current_layer[j]
            next_point = next_layer[j]

            all_lines.append((current_point, next_point))

            prev_pos_index = (j - 1) % len(current_layer)  # Wrap around to the last node
            next_point_prev_pos = next_layer[prev_pos_index]
            all_lines.append((current_point, next_point_prev_pos))

            next_pos_index = (j+1)%len(current_layer)
            next_pos = next_layer[next_pos_index]
            all_lines.append((current_point, next_pos))

create_lines2(all_nodes, [0,1])
create_lines2(all_nodes, [1,2])
create_lines2(all_nodes, [2,3])
create_lines2(all_nodes, [3,4])
create_lines2(all_nodes, [4,5])


unique_lines = set()

# Iterate through the lines and add them to the set
for line in all_lines:
    a, b = line
    if (a, b) not in unique_lines and (b, a) not in unique_lines:
        unique_lines.add((a, b))

# Convert the set back to a list if needed
unique_lines_list = list(unique_lines)
sorted_list = sorted(unique_lines_list, key=lambda line: line[0])


# Specify the file name
lines_file = 'lines.txt'

# Open the file for writing
with open(lines_file, 'w') as file:
    for line in sorted_list:
        point1, point2 = line
        x1, y1, z1 = point1
        x2, y2, z2 = point2
        file.write(f'({x1:.2f}, {y1:.2f}, {z1:.2f}) - ({x2:.2f}, {y2:.2f}, {z2:.2f})\n')

#print(f'Lines saved to {output_file}')


#Not in use
def midpoints(lines):
    midpoints = []
    for line in lines:
        point1, point2 = line
        x1, y1, z1 = point1
        x2, y2, z2 = point2
    
        # Calculate the midpoint coordinates
        mid_x = (x1 + x2) / 2.0
        mid_y = (y1 + y2) / 2.0
        mid_z = (z1 + z2) / 2.0
    
        midpoint = (mid_x, mid_y, mid_z)
        midpoints.append(midpoint)
    return midpoints



'''
def horizontal_lines(nodes,layer):
    my_list = nodes[layer]
    for i in range(len(my_list)):
        point1 = my_list[i]
        point2 = my_list[(i + 1) % len(my_list)]  # Circular pairing
        all_lines.append((point1, point2))
        

def vertical_lines(nodes, layer1, layer2):
    for i in range(len(nodes[layer1])):
        point1 = nodes[layer1][i]
        point2 = nodes[layer2][i]
        all_lines.append((point1,point2))
        

def create_lines(nodes, layer_indices):
    for i in range(len(layer_indices) - 1):  # Iterate through specified layers
        current_layer = nodes[layer_indices[i]]
        next_layer = nodes[layer_indices[i + 1]]

        for j in range(len(current_layer)):
            current_point = current_layer[j]
            next_point = next_layer[j]

            all_lines.append((current_point, next_point))
            
            prev_pos_index = (j - 1) % len(current_layer)  # Wrap around to the last node
            next_point_prev_pos = next_layer[prev_pos_index]
            all_lines.append((current_point, next_point_prev_pos))

'''