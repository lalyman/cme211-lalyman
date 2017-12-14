import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import scipy.sparse
import scipy.sparse.linalg

class Truss:
    """Returns a data table of forces from each beam.

    Specifically, this class creates a system of
    linear equations by balancing the x and y forces
    acting on each joint, where the forces provided
    by the beams are variables. Then the class solves
    this system of equations to determine each beam
    force and outputs the results in a data table.
    If the system cannot be solved (i.e. there is
    a singular matrix or non-square matrix) the class
    returns an error with an explanation.

    Also, the code creates a figure of the truss system;
    that is, the figure has a picture of the joints and 
    the beams connecting them, regardless of whether the
    corresponding matrix equation has a solution.
    The user can input an (optional) file path, and the code
    saves the figure in that file path location.
    """
    def __init__(self,joints_file,beams_file,plot_file):
        """Initializes the class and its methods.

        Gives the Truss class attributes 'joints_file'
        and 'beams_file' from the user-inputted data
        files. Also gives the Truss class a 'plot_file'
        attribute, which is the plot's name/where it
        should be saved if the user specifies a file path;
        if no file path is specified, the plot is is named
        "figure" by default and saved in the working
        directory. otherwise. Then the class methods
        are initialized."""
        self.joints_file = str(joints_file)
        self.beams_file = str(beams_file)
        self.plot_file = str(plot_file)
        self.get_beam_data(self.beams_file)
        self.get_joints_data(self.joints_file)
        self.PlotGeometry()
        self.get_beam_length()
        self.check_if_square()
        self.get_beam_eqn()
        self.make_b_vec()
        self.solve_linear_system()

    def get_beam_data(self,beams_file):
        """Gets information from beam data file.

        Creates the 'beam_joints' attribute,
        which is a dictionary where the keys
        are beams and the values are the joint
        indices of the two joints per each beam."""
        self.beam_joints = {}
        with open(beams_file,'r') as f:
            # Ignore the first line (which
            # is just the table heading)
            header = f.readline()
            # Iterate through the file with
            # beam data
            for line in f:
                line = line.split()
                beam = int(line[0])
                # Remove the two joint indices
                # for the beam
                Ja = int(line[1])
                Jb = int(line[2])
                self.beam_joints[beam] = [Ja,Jb]


    def get_joints_data(self,joints_file):
        """Gets information from joints data file.

        In particular, this method adds three class
        attributes that are all dictionaries:
        joints_with_attached_beams, joint_position,
        and joint_forces. All of the dictionaries
        have the joints (i.e. joint indices) as keys.
        Dictionary 'joint_position' has values that
        are the (x,y) coordinates for each joint. Then
        'joint_forces' has values that are lists
        including: the x-forces on each joint, the y-forces
        on each joint, and a boolean for whether each joint
        is attached to a wall. Finally, 'joints_with_attached
        _beams' has values that are the beam indices of the
        beams attached to each joint."""
        self.joints_with_attached_beams = {}
        self.joint_position = {}
        self.joint_forces = {}
        with open(joints_file,'r') as f:
            # Skip the header line
            header = f.readline()
            # Iterate through the file
            for line in f:
                line = line.split()
                joint = int(line[0])
                # Collect the data
                x = float(line[1])
                y= float(line[2])
                Fx = float(line[3])
                Fy = float(line[4])
                is_attached = int(line[5])
                # Add data to dictionaries
                self.joint_position[joint] = [x,y]
                self.joint_forces[joint] = [Fx,Fy,is_attached]
                # Calls helper function to add data to
                #'joints_with_beams attached' dictionary
                self.add_to_joints_with_attached_beams_dict(joint)
        f.close()

    def add_to_joints_with_attached_beams_dict(self,joint):
        """Collects beam indices of beams attached to inputted joint.

        This is a helper function used by the 'get_joints_data' method."""
        # Iterates through existing beam_joints dictionary
        for beam in self.beam_joints.keys():
            # Checks if joint is attached to the beam
            if joint in self.beam_joints[beam]:
                # Check if the inputted joint is already a key in
                # the dictionary. If not, creates a new dictionary
                # entry
                if joint not in self.joints_with_attached_beams.keys():
                    # Add new dictionary entry
                    self.joints_with_attached_beams[joint] = [beam]
                # If the joint is already a key, adds the beam to the
                # list of beams attached to that joint 
                else:
                    self.joints_with_attached_beams[joint].append(beam)

    def get_beam_length(self):
        """Computes the length of a beam and some trig. information.

        Calls its helper method 'get_trig_info' to
        populate the dictionary attributes 'beam_length'
         and 'beam_trig'. For 'beam_trig' 
        """
        self.beam_lengths = {}
        self.beam_trig = {}
        for beam in self.beam_joints.keys():
            joints = self.beam_joints[beam]
            # Select out each joint
            Ja = joints[0]
            Jb = joints[1]
            pos_a = self.joint_position[Ja]
            pos_b = self.joint_position[Jb]
            # Calls helper function for computations
            [dist,cos_angle,sin_angle] = self.get_trig_info(pos_a,pos_b)
            self.beam_lengths[beam] = dist
            # Trig values are relative to Ja (first joint)
            self.beam_trig[beam] = [cos_angle,sin_angle]

    def get_trig_info(self,point1,point2):
        """Computes trigonometric info for two points. 

        Specifically, computes the Euclidean distance
        between the inputted points and the cosine/
        sine of the angle formed by viewing
        the points as vectors. This is a helper
        method called by the method 'get_beam_length'.
        """
        deltax = point2[0] - point1[0]
        deltay = point2[1] - point1[1]
        # Get Euclidean distance
        dist = math.sqrt(deltax**2 + deltay**2)
        # Get angle information (relative to first point)
        cos_angle = deltax/dist
        sin_angle = deltay/dist
        return [dist,cos_angle,sin_angle]

    def check_if_square(self):
        """Checks if the matrix of data is square.

        Specifically, the method throws an error if
        the matrix is rectangular i.e. the number of
        equations (m) does not equal the number of
        variables (n)."""

        # Two equations per joint, and the number of joints is
        # the number of keys in the dictionary 'joint_position'
        self.number_of_eqns = len(self.joint_position.keys())*2
        self.number_of_beams = len(self.beam_trig.keys())
        number_attached_joints = 0
        not_square_error = "Linear system is not square. Method of " + \
            "joints isn't suitable."
        # Iterate through the joints
        for val in self.joint_forces.values():
            is_attached = val[2]
            # Get total number of attached joints
            number_attached_joints += is_attached
        # There is one variable per beam. Also, there are two
        # additional variables per attached joint representing
        # the force of attachment in the x and y directions.
        self.number_of_vars = self.number_of_beams + 2*number_attached_joints
        # Check if matrix is not square
        if self.number_of_vars != self.number_of_eqns:
            raise RuntimeError(not_square_error)

    def update_CSR_arrays(self,joint,beams_attached,val_array,col_array,
            row_array,number_of_attached,is_attached,direc_bool):
        for beam in beams_attached:
            # Direction is a boolean for which 0 indicates
            # the 'x-direction' and 1 indicates the 'y-direction'
            joints_of_beam = self.beam_joints[beam]
            # Check if joint is the first joint of the beam
            if joint == joints_of_beam[0]:
                # If so, get relevant angle info
                beam_angles = self.beam_trig[beam]
           # Otherwise, joint must be the second joint of the beam
            else:
                # Multiply angles by -1 since trig. is relative
                # to the first joint
                beam_angles = [-1*i for i in self.beam_trig[beam]]
            # Build up arrays for CSR format
            if beam_angles[direc_bool] != 0:
                val_array.append(beam_angles[direc_bool])
                # Indexing starts at 0
                col_array.append(beam - 1)
            # Add contribution from attachment force (if one exists)
            if is_attached != 0:
                # Always '1' added
                val_array.append(1)
                col_array.append(self.number_of_beams+number_of_attached+
                    direc_bool)
        # Equation done, switching to a new row of A
        row_array.append(len(val_array))
        return [val_array,col_array,row_array]

    def get_beam_eqn(self):

        val_array = []
        col_array = []
        # CSR row array always starts with a zero
        row_array = [0]
        number_of_attached = 0
        # Iterate through the joints
        for joint in self.joints_with_attached_beams.keys():
            beams_attached = self.joints_with_attached_beams[joint]
            is_attached = self.joint_forces[joint][2]
            # For the joint, add equation for balancing x-direction    
            [val_array,col_array,row_array] = self.update_CSR_arrays(
                joint,
                beams_attached,
                val_array,
                col_array,
                row_array,
                number_of_attached,
                is_attached,
                0)     
             # For the joint, add equation for balancing the y-direction
            [val_array,col_array,row_array] = self.update_CSR_arrays(
                joint,
                beams_attached,
                val_array,
                col_array,
                row_array,
                number_of_attached,
                is_attached,
                1)
            # Update the number of equations arising from
            # joints being attached to walls (2 equations, one
            # for x-direction and one for y-direction, per joint
            # attached to a wall)
            if is_attached != 0:
                number_of_attached +=2

        A = scipy.sparse.csr_matrix((val_array,col_array,row_array))
        self.check_for_singular_matrix(A)
        self.A = A

    def check_for_singular_matrix(self,A):
        """Throws error if linear system from data files has singular matrix.

        Inputted matrix A is in CSR (compressed sparse row) format."""
        error_message = "Matrix is singular, so linear system " + \
            "cannot be solved."
        # Converts CSR to CSC to remove extraneous
        # Python warning message when using 'linalg.inv'
        A = scipy.sparse.csc_matrix(A)
        try:
            inv = scipy.sparse.linalg.inv(A)
        except RuntimeError as err:
            raise RuntimeError(error_message)

    def make_b_vec(self):
        """Creates b vector for solving Ax = b.

        The entries of b are the external x/y forces
        acting on each joint that aren't from beams
        i.e. forces from being attached to a wall.
        An entry of 0 indicates that the corresponding
        joint is either: (1) not attached, or (2) attached
        in a way that doesn't produce force in that
        direction. Vector b is outputed in CSR
        format."""
        # Build arrays for CSR format
        col_array = []
        val_array = []
        row_array = [0]
        for forces in self.joint_forces.values():
            # Select out forces in each direction
            force_x = forces[0]
            force_y = forces[1]
            # These checks for 0 are purely for formatting,
            # so values of 0 are not assigned '+' or '-' signs
            if force_x != 0:
                val_array.append(force_x)
                # Vector b is a single column, so '0' is 
                # always the column index
                col_array.append(0)
            row_array.append(len(val_array))
            if force_y != 0:
                val_array.append(force_y)
                # Vec. b is a single column, so '0' is column index
                col_array.append(0)
            row_array.append(len(val_array))
        # Use arrays to get b in CSR format
        b = scipy.sparse.csr_matrix((val_array,col_array,row_array), \
            shape=(self.number_of_eqns, 1))
        return b 

    def solve_linear_system(self):
        A = self.A
        b = self.make_b_vec()
        soln = scipy.sparse.linalg.spsolve(A, b)
        force_soln = {}
        for i in range(self.number_of_beams):
            # Remove signed zeros
            if soln[i] == 0: soln[i] = abs(soln[i])
            # Adjust beam number up by one index
            force_soln[i+1] = soln[i]
        self.force_soln = force_soln

    def PlotGeometry(self):
        """Creates (optional) figure of the beam/truss system."""
        file = self.plot_file
        plt.figure(1)
        xvals = []
        yvals = []
        for beam in self.beam_joints.keys():
            Ja = self.beam_joints[beam][0]
            Jb = self.beam_joints[beam][1]
            Ja_pos = self.joint_position[Ja]
            Jb_pos = self.joint_position[Jb]
            xvals += [Ja_pos[0], Jb_pos[0]]
            yvals +=[Ja_pos[1], Jb_pos[1]]
            plt.plot(xvals,yvals,'b-')
        # Improve formatting for plot
        xmin = min(xvals)
        xmax = max(xvals)
        ymin = min(yvals)
        ymax = max(yvals)
        padding = .05 *(xmax - xmin)
        plt.tight_layout()
        plt.show(block=False)
        # Saves the file. All file handling issues are checked
        # ahead of time in 'main.py'
        plt.savefig(file)
   
    def make_table_header(self):
        """Creates a table header for outputted data table.

        This is a helper function used by '__repr__' below."""
        string = " Beam       Force\n"
        string += "-----------------\n"
        return string

    def __repr__(self):
        """Returns final data table to console."""

        # Call helper function to make the table header
        table_string = self.make_table_header()
        colsize = 9
        # Populate the data table and format it
        for beam in self.force_soln:
            line = '    %i   %%' % beam
            line += '%i' % colsize + '.3f' 
            line = line % self.force_soln[beam] + '\n'
            table_string += line 
        # Output to the console
        return table_string
