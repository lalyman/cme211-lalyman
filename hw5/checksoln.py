import numpy as np
import os.path
import sys

# Provide usage message
if len(sys.argv) != 3:
    print('Usage:')
    print('  python3 %s <maze file> <solution file>'% sys.argv[0])
    sys.exit(0)

maze_file = str(sys.argv[1])
soln_file = str(sys.argv[2])

# Error messages for inputted files
maze_file_not_found = "ERROR: The maze file '%s' was not found." \
% maze_file
soln_file_not_found = "ERROR: The solution file '%s' was not found."\
% soln_file
# Error messages for incorrect solution
bad_entrance = "ERROR: The solution file " + \
"'%s' does not start at the maze entrance." % soln_file
not_done = "ERROR: The solution file " + \
"'%s' is incorrect, because it never exits the maze." % soln_file
illegal_move = "ERROR: The solution file moves more than one " + \
"position in one step. You can only move 1 unit north, south, " + \
"east, or west per move."

# Check if user-inputted files are there

if os.path.isfile(soln_file) is False:
    print(soln_file_not_found)
    sys.exit(0)

# Read the maze file and store the data
with open(maze_file,'r') as f:
      # Use first line to get number of rows
      # and columns
    dimensions = f.readline()
    dimensions = dimensions.split()
    num_rows = int(dimensions[0])
    num_cols = int(dimensions[1])
    # Initialize 2D numpy array 'maze_data' to all zeros
    # Will add 1s to indicate presence of a wall
    # (Instructions said to use numpty array)
    maze_data = np.zeros((num_rows,num_cols))
    # Iterate through the maze file and
    # store wall positions in 'maze_data'
    for line in f:
        line = line.split()
        row = int(line[0])
        col = int(line[1])
        # Set to 1 to indicate that there is a wall
        # in this position
        maze_data[row][col] = 1
f.close()
# Initialize array for storing data from user-inputted solution file
soln_data = []
entrance = -1
# Get true starting position in maze
for j in range(num_cols):
    if maze_data[0][j] != 1:
        entrance = j
# Now read through the solution file to check if 
# solution is valid
with open(soln_file,'r') as f:
    # Get the starting point
    first_pos = f.readline()
    first_pos = first_pos.split()
    first_row = int(first_pos[0])
    first_col = int(first_pos[1])
    soln_data.append([first_row,first_col])
    # Check if maze was properly entered on the first row    
    if first_row != 0 or first_col != entrance:
        print(bad_entrance)
        sys.exit(0)
    # Iterate through the rest of the positions in the
    # solution file
    for line in f:
        line = line.split()
        row = int(line[0])
        col = int(line[1])
        soln_data.append([row,col])
        # Check if you went out of the bounds of the maze
        if row >= num_rows or col >= num_cols or \
        row < 0 or col < 0:
            out_of_bounds = "ERROR: The solution file has an entry "+ \
            "(%s,%s) that is out of the maze boundaries." %(row, col)
            print(out_of_bounds)
            sys.exit(0)
        # Check if you hit a wall
        if maze_data[row][col] == 1:
            hit_a_wall = "ERROR: The solution file has an entry " + \
            "(%s,%s) that hits a wall!" % (row, col)
            print(hit_a_wall)
            sys.exit(0)
f.close()

last_pos = soln_data[-1]
# Check if the exit of the maze was reached on the last row.
if last_pos[0] != (num_rows - 1):
    print(not_done)
    sys.exit(0)
# Ensure each position change was valid (i.e. you moved one 
# position at a time) by iterating through 'soln_data'
for i in range(len(soln_data)-1):
    [row_i,col_i] = soln_data[i]
    [row_next,col_next] = soln_data[i+1]
    delta_x = abs(row_next-row_i)
    delta_y = abs(col_next-col_i)
    if delta_x != 1 and delta_y != 1:
        print(illegal_move)
        sys.exit(0)
# If no error was thrown by now, the solution must be valid
print('Solution is valid!')
