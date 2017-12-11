CME 211: Homework 3
Laura Lyman

---------------------------------------
DESCRIPTION OF ASSIGNMENT:
In this assignment, we create an Airfoil class to process and display
data about one or more airfoils, which are 2D cross sections of aircraft
wings. Specifically, given wind flow at some angle of attack (alpha) with 
respect to an airfoil's leading edge, we are trying to:
-- compute the associated lift coefficient c_l, 
-- locate any stagnation points, which are points (x,y) on this 2D cross 
section at which the flow velocity is zero, and
-- determine the corresponding pressure coefficient C_p for each
stagnation point found.
All of this information should be displayed in a table in the console for
each angle alpha for which we have data.

The airfoil geometry is defined by a series of points, where
each pair of consecutive points can be connected by a straight line to 
form 'panels' of the airfoil. In a given directory (e.g. 'naca0012' or
'naca2412'), pressure coefficient data is available for each of these
panels for various angles of attack (alpha), where each 'alpha file' 
contains data for all of the panels for that particular angle. In addition,
a directory should include an x-y data file (i.e. 'xy.dat') that contains 
the series of points defining the airfoil panels. With this information, 
the pressure coefficient distribution can then be integrated across the 
airfoil surface to determine the lift coefficient c_l. The stagnation 
point(s) occur where the pressure coefficient C_p is 1.0, so they can be
found by simply reading through a given alpha file; if no point has a 
pressure coefficient of exactly 1.0, we simply return the panel whose C_p
is closest to 1.0 and an (x,y) coordinate that is the average of the points
defining that panel.

---------------------------------------
OOP CONCEPTS:
In object oriented programming, code structure is based around the 
concept of 'objects', which are instances of certain classes. In 
particular, classes define the properties and behavior of the 
objects they represent. Objects have two kinds of attributes 
(accessed by the '.'syntax): data attributes (instance variables) 
and function attributes (methods). For example, in my Airfoil class, 
there is the 'dir' or directory attribute, accessed by '.dir' 
(specifically by self.dir within the class definition). Similarly, 
there is a 'check_valid_directory()' function attribute of Airfoil that 
is used to determine whether the user inputted a valid directory for the 
class to use i.e. if 'self.dir' refers to a sensible file path.

In the Airfoil class, the instance variables are:

self.dir 			(user inputted file path of directory)
self.alpha_files	(dictionary where the values are the
					paths of the 'alpha<angle of attack>.dat'
					files and the keys are the corresponding
					angle of attack values e.g. -3.0, 0.0,
					etc.)
self.xy_data_file		(file path of 'xy.dat')
self.chord_len			(length of chord across airfoil)
self.x_y_vals			(dictionary whose values are (x,y)
				coordinates from 'xy.dat' and keys
				are the unique corresponding line
				numbers for each coordinate in the
				file, where numbering starting from 0
				on the first line of data after the
				header line)
self.alpha_stag_dict		(dictionary whose keys are alpha values
				and whose values are lists. A list's
				first entry is the x-y coordinates of a
				stagnation point for that alpha key
				and its second entry is the pressure 
				coefficient C_p associated with the
				corresponding stagnation point, repeating
				in this pattern if an angle of attack has 
				more than one stagnation point)
self.deltas 			(dictionary whose keys are line indices
				in 'xy.dat' corresponding to certain
				pairs of panels and whose values are
				each [delta_x,delta_y], which is the 
				distance differences in the x/y direct-
				ions between that pair of panels)
self.line_label			(header label on the first line of 'xy.
				dat' (e.g. 'NACA 0012'), which is used in
				the console output at the end)
						
and the 12 function attributes are:

-- self.check_valid_directory()
-- self.check_for_needed_files()
-- self.get_alpha_val_from_local_file(local_file)
-- self.compute_xy_vals()
-- self.check_data_format_xy_file(data)
-- self.get_table_data(file,alpha)
-- self.get_stag_point(Cp_val,stag_cp_val,i,x,y)
-- self.build_divider(number)
-- self.add_two_lists_elementwise (self,list1,list2)
-- self.convert_deg_to_rad(angle)
-- self.compute_deltas()
-- self.make_console_output_header()

along with __init(self,directory) and __repr__(self). The document strings 
for all of these methods provide more specifics.

Three central concepts of OOP are:

ABSTRACTION: representing data/computations in a familiar form (making
code more human-readable). The Airfoil class uses abstraction by having
helper functions (with very human-readable/understandable names like
'convert_deg_to_rad()' or 'compute_xy_vals()') do computations on
the side that would otherwise interrupt the code flow or distract from the 
code's general structure. In addition, data attributes of the Airfoil class 
are given intuitive names to make the code more understandable; for example,
the data attribute 'self.chord_len' represents the airfoil's chord length.

ENCAPSULATION: hiding the details of data structures and algorithms 
(internal code). That is, encapsulation is a way to prevent the user from 
changing the state of an object in a way that is not desirable. One way 
to implement encapsulation is to make methods more "private" by telling 
Python to obscure/mangle attribute names, which is done by placing '__' 
before an attribute definition. For example, my Airfoil class protects
the data initialization and final console output with '__init(self,dir)'
and '__repr__(self)'. 

Another way to demonstrate encapsulation is to return references/copies of
class attributes that *can't* be changed by a user rather than the (mutable)
attributes themselves. For example, instead of letting a user access/change 
the Airfoil's chord length directly, we would provide a method 

def get_chord_len(self):
	return copy.deepcopy(self.chord_len)

so that the user would interact with an independent chord length that would
not affect the airfoil's actual chord length. However, since this assignment 
is concerned more with displaying a data table than with having a user 
interact with instances of the Airfoil class, these types of protective
'get_some_attribute()' methods were not used.

DECOMPOSITION: breaking large problems into smaller, more manageable
problems. The Airfoil class uses decomposition with its large number
of methods/helper functions.

---------------------------------------
ERROR HANDLING:

RuntimeError error exceptions are used to check:

-- if the directory (file name) given by the user is
a real directory
-- if the user-inputted directory has at least one file
of the form 'alpha<angle of attack>.dat' i.e. if
at least one file starts with 'alpha' and is 
followed by an <angle of attack> that is a float
-- if the user-inputted directory has an 'xy.dat' file
-- if the data in 'xy.dat' has the correct form (i.e.
after the header, all entries are floats) so the data
can actually be read
-- if the data in any of the alpha files has the correct
form (i.e. beyond the header row, all entries are floats).

In addition, a useful 'Usage' message is printed when the user
inputs an insufficient number of arguments (i.e. doesn't provide a 
directory).

---------------------------------------
RESULTS:
The method was then run on the 'naca0012' directory to produce:

Test case: NACA 0012

   alpha      cl           stagnation pt
   -----   -------   --------------------------
   -3.00   -0.3622   ( 0.0030,  0.0094)  0.9906
    0.00    0.0000   ( 0.0000,  0.0000)  0.9944
    3.00    0.3622   ( 0.0030, -0.0094)  0.9906
    6.00    0.7235   ( 0.0099, -0.0170)  0.9967
    9.00    1.0827   ( 0.0219, -0.0246)  0.9977


matching the output on the bottom of page 4 in the
HW3 PDF. Similarly, on the 'naca2412' folder we have
the output:

Test case: NACA 2412

   alpha      cl           stagnation pt
   -----   -------   --------------------------
   -3.00   -0.1072   ( 0.0021,  0.0081)  0.9926
    0.00    0.2554   ( 0.0000,  0.0005)  0.9936
    3.00    0.6172   ( 0.0027, -0.0088)  0.9917
    6.00    0.9774   ( 0.0097, -0.0158)  0.9966

As a sanity check, note that each stagnation point
occurs when the associated pressure coefficient C_p 
(last entry in each row) is nearly 1.0, as expected.
