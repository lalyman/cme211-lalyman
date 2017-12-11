import glob
import math
import os.path

class Airfoil:
    """Returns a data table.

    Specifially, the method returns a table of
    the lift coefficient (cl),stagnation point,
    and pressure coefficient for various angles
    of attack (alpha) given alpha data files in
    the user-inputted directory."""
    def __init__(self,directory):
        """Initializes the class and its methods.

        Gives the class a 'directory' attribute
        from the user-inputted file path to a
        directory."""
        self.dir = str(directory)
        self.check_valid_directory()
        self.check_for_needed_files()
        self.compute_x_y_vals()
        self.compute_deltas()

    def check_valid_directory(self):
        """Checks if inputted directory exists.

        Throws an error if the directory is not 
        found. Note that the empty string is
        *not* considered a valid directory.
        Then removes the final '/' character
        if it's present at the end of the
        directory file path so that formatting
        is consistent."""
        directory_not_found = "Directory not found."

        if os.path.isdir(self.dir) is False:
            raise RuntimeError(directory_not_found)

        if self.dir.endswith('/'):
            self.dir = self.dir[:-1]

    def check_for_needed_files(self):
        """Checks if directory has necessary files.

        Specifically, the method determines whether
        the inputted directory contains an 'xy.dat'
        file and at least one file of the form
        'alpha<angle of attack>.dat' and throws an
        error if these files are not present. 
        The alpha file(s) and xy data file are then
        stored as dictionary attributes for later use.""" 
        no_alpha_file = "Directory does not have a " +\
                "file of the form 'alpha<angle of attack>.dat'"
        no_xy_file = "Directory does not have an " +\
                "'xy.dat' file."
        has_alpha_file = False
        has_xy_file = False
        self.alpha_files = {}
        # Uses os.path to iterate through all the
        # items in the directory.
        for path in glob.glob(self.dir + '/*'):
            # Gets the file name relative to being
            # in the directory.
            local_file = os.path.split(path)[1]
            # Checks for alpha file
            if local_file.startswith("alpha"):
                has_alpha_file = True
                alpha = self.get_alpha_val_from_local_file(local_file)
                # Stores the file path as a value
                # in the 'alpha_files' dictionary
                # with key equal to the angle of
                # attack (alpha)
                self.alpha_files[alpha] = path
            # Checks for 'xy.dat' file and adds it
            # as a class attribute
            if local_file.startswith("xy"):
                has_xy_file = True
                self.xy_data_file = path
        if has_alpha_file is False:
            raise RuntimeError(no_alpha_file)
        if has_xy_file is False:
            raise RuntimeError(no_xy_file)

    def get_alpha_val_from_local_file(self,local_file):
        """Helper function to get <angle of attack>.

        Inputted 'local_file' has the form 'alpha<angle 
        of attack>.ext' where '.ext' is some extension.
        For example, inputting 'alpha-3.0.dat' returns
        the float -3.0. Function is used by the 
        check_for_needed_files() method above."""
        prefix = "alpha"
        bad_file_name = "Alpha file not present or file(s) are " \
                "named incorrectly. You need at least " + \
                "one alpha file with the form " + \
                "'alpha<angle of attack>.dat' " + \
                "where <angle of attack>  is a float."

        file_without_ext= os.path.splitext(local_file)[0]
        alpha = file_without_ext[len(prefix):]
        # Extra error checking: ensures the angle of
        # attack is a float i.e. if the file is properly
        # named.
        try:
            alpha = float(alpha)
        except ValueError:
            raise RuntimeError(bad_file_name)
        return float(alpha)

    def compute_x_y_vals(self):
        """Returns (x,y) coordinates of each panel.

        Reads 'xy.dat' to get x-y values for each
        panel and checks if the data is in the 
        correct format (i.e. values can be converted
        to floats). Also calculates the airfoil 
        chord length."""
        x_y_vals = {}
        # Bogus initial values
        x_min = 2
        x_max =- 2
        
        with open(self.xy_data_file,'r') as f:
            # Take out the fist (text) line so
            # it can be outputted to the console
            # at the end. The coordinate data
            # starts on the second line.
            self.line_label = f.readline()
            index = 0
            for line in f:
                data = line.split()
                # Checks if coordinate data is
                # formatted correctly.
                self.check_data_format_xy_file(data)
                # Puts (x,y) values in dictionary 'x_y_vals'.
                x_y_vals[index] = data
                x_y_vals[index][0] = float(x_y_vals[index][0])
                x_y_vals[index][1] = float(x_y_vals[index][1])
                # Calculates x_min and x_max to determine
                # chord length.
                if x_y_vals[index][0] >= x_max: x_max = x_y_vals[index][0]
                if x_y_vals[index][0] <= x_min: x_min = x_y_vals[index][0]
                index += 1
        f.close()
        self.chord_len = abs(x_max - x_min)
        self.x_y_vals = x_y_vals

    def check_data_format_xy_file(self,data):
        """Check if 'xy.dat' has C_p data in correct format.

        The variable 'data' is a list of the entries in a
        particular line of the xy file. The method checks if
        the data values are indeed floats or can be 
        converted to floats without issue."""
        data_format_error = "Data in 'xy.dat' not in correct " + \
                            "format."
        for val in data:
            try:
                val = float(val)
            except ValueError:
                raise RuntimeError(data_format_error)       

    def get_table_data(self,file,alpha):
        """Yields stag. pts with lift/pressure coefficients.

        The method computes these values (cl,cp, and the 
        stagnation point) for each angle of attack alpha 
        given by each alpha data file."""
        with open(file,'r') as f:
            i = 0
            c_x = 0
            c_y = 0
            # Dictionary with alpha values as keys
            # (e.g. -3.0, 0.0, etc) and two element
            # lists as values. The first entry of
            # each list value is the (x,y) coordinates
            # of the stagnation point and the second 
            # entry is the pressure coefficient C_p 
            # of the stagnation point.
            self.alpha_stag_dict = {}
            # Skip first line (header)
            header_line = next(f)
            # Bogus values
            stag_cp_val = -5
            [x,y] = [-2,-2]
            format_error = "Data in 'alpha<angle of attack>.dat' "+ \
                "not in correct format."
            for line in f:
                # Check if data in alpha file has
                # correct format (i.e. is a float)
                try: Cp_val = float(line.split()[0])
                except ValueError: raise RuntimeError(format_error)
                Cp_val = float(line.split()[0])
                # Compute stagnation point
                [stag_pt,stag_cp_val]= self.get_stag_point(Cp_val,\
                stag_cp_val,i,x,y)
                [x,y] = stag_pt
                # Compute differences in x and y
                [deltax, deltay] = self.deltas[i]
                c_x += -1*Cp_val*deltay / (self.chord_len)
                c_y += Cp_val*deltax / (self.chord_len)
                i += 1
        f.close()
        # Convert angle of attack to radians!
        angle = self.convert_deg_to_rad(alpha)
        self.alpha_stag_dict[alpha] = [(x,y), stag_cp_val]
        # Compute lift coefficient
        cl = c_y*math.cos(angle) - c_x*math.sin(angle)
        return cl

    ##### Begin section of (smaller) helper functions #####

    def get_stag_point(self,Cp_val,stag_cp_val,i,x,y):
        """Calculates/adjusts stagnation pt and C_p.

        At a point where the pressure coefficient C_p
        is sufficiently close to 1.0, the method returns
        the panel value closest to 1.0 and an (x,y)
        coordinate that is the average of the points
        defining that particular panel. Used in method
        get_table_data() for computations as alpha
        files are read."""

        # Get C_p value closest to 1
        if abs(Cp_val - 1) <= abs(stag_cp_val - 1):
            stag_cp_val = Cp_val
            # Take average of points
            pt = self.add_two_lists_elementwise( \
                self.x_y_vals[i],self.x_y_vals[i+1])
            x = 1/2*pt[0]
            y = 1/2*pt[1]
        return [[x,y],stag_cp_val]


    def build_divider(self,number):
        """Outputs a string of '---' characters.

        The number of '-' characters in this line
        can be specified. Used for formatting the
        console output."""
        s = ""
        for i in range(number):
            s += "-"
        return s

    def add_two_lists_elementwise(self,list1,list2):
        """Adds two lists elementwise & returns new list.

        Used in method get_stag_point to average the 
        coordinates of two panels."""
        newlist = []
        for i in range(len(list1)):
            newentry = list1[i] + list2[i]
            newlist += [newentry]
        return newlist
    
    def convert_deg_to_rad(self, angle):
        """Converts attack angle from deg. to radians."""
        rad = angle * 2 * math.pi / 360
        return rad

    def compute_deltas(self):
        """Computes x and y distances btwn. pairs of panels."""
        x_y_vals = self.x_y_vals
        # Stores delta(x) and delta(y) for each pair in
        # a dictionary, where the keys are indices
        # referring to certain pairs given the order 
        # in which the panels are traversed.
        delta_dict = {}
        # For given data, there are 160 (x,y) values
        number_of_vals = len(x_y_vals.keys())
        for i in range(0,number_of_vals-1):
            deltas1 = x_y_vals[i]
            deltas2 = x_y_vals[i+1]
            deltax = deltas2[0] - deltas1[0]
            deltay = deltas2[1] - deltas1[1]
            delta_dict[i] = [deltax,deltay]
        self.deltas = delta_dict

    def make_console_output_header(self):
        """Creates table header for output to console."""
        divider = "-------- "
        string = "Test case: " + self.line_label
        # Make columns
        string += '\n' + "%8s" % "alpha" + "%8s" % "cl" \
                + "%24s" % "stagnation pt"
        string += '\n' + "%8s" % self.build_divider(5) + \
                "%10s" % self.build_divider(7) + \
                "%29s" % self.build_divider(26)
        return string

    ##### Output data table to console #####

    def __repr__(self):
        """Populates the data table."""
        # Make header of outputted table
        table_string = self.make_console_output_header()
        for key in sorted(self.alpha_files):
            path = self.alpha_files[key]
            cl = self.get_table_data(path,key)
            # Get stagnation point and its pressure
            # coefficient C_l
            stag_x = "%7.4f" % self.alpha_stag_dict[key][0][0]
            stag_y = "%7.4f" % self.alpha_stag_dict[key][0][1]
            stag_val = "%8.4f" % self.alpha_stag_dict[key][1]
            stagnation_pt = "%21s" % ('(' + stag_x + ', ' + \
                stag_y + ')')
            # Put data into table
            table_string += '\n' + "%8.2f" % key + "%10.4f" % cl \
                + stagnation_pt + stag_val
        return table_string
