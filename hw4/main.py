import os.path
import sys
import truss

### BEGIN HELPER FUNCTIONS ###
def bad_format(file):
    format_error = "ERROR: Make sure \'{0}\'".format(file)
    format_error += " does not start or end with a \'/\' character."
    if file[0] == '/' or file[-1] == '/':
        print(format_error)
        return True

def file_not_found(file):
    error = "ERROR: File '%s' not found." % file
    if os.path.isfile(file) is False:
        print(error)
        return True

def directory_not_found(directory):
    error = "ERROR: Directory '%s' not found." % directory
    if os.path.isdir(directory) is False:
        print(error)
        return True
### END HELPER FUNCTIONS ###

# Prints usage message
if len(sys.argv) < 3:
    print('Usage:')
    print('  python3 %s [joints file] [beams file] \
[optional plot output file]' % sys.argv[0])
    sys.exit(0)

joints_file = sys.argv[1]
beams_file = sys.argv[2]
# Checks for formatting issues in user inputted files
if bad_format(joints_file) or bad_format(beams_file):
    sys.exit(0)
# Checks if user inputted files exist
if file_not_found(joints_file) or file_not_found(beams_file):
    sys.exit(0)

# Handles when user does not submit the (optional) plot file name
# Just sets the name of the plot to be "figure.pdf" and saves it in
# the working directory
if len(sys.argv) == 3:
    try:
        t = truss.Truss(joints_file,beams_file,"figure")
    except RuntimeError as e:
        print('ERROR: {}'.format(e))
        sys.exit(2)
    print(t)

# Handles when user does include the output plot file name
# Picks third argument as plot file, even if more than 3 arguments
# are inputted (for some reason)
if len(sys.argv) > 3:
    plot_file = sys.argv[3]
    # Check for formatting issues on user inputted plot file
    if bad_format(plot_file):
        sys.exit(0)
    [plot_filepath,plot_filename] = os.path.split(plot_file)
    # If user attached a directory to their plot file name, check if
    # that directory exists
    # Otherwise the plot file will just be saved in the current directory
    if len(plot_filepath) != 0 and directory_not_found(plot_filepath):
        sys.exit(0)
    
    try:
        t = truss.Truss(joints_file,beams_file,plot_file)
    except RuntimeError as e:
        print('ERROR: {}'.format(e))
        sys.exit(2)
    print(t)
