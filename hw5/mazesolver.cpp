#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::get;
using std::make_tuple;
using std::string;
using std::tuple;
// Define constants to size static array
// and set available storage
#define num_rows_avail 500
#define num_cols_avail 500

// Problem assumptions:
//
// -- each maze has a solution
// -- the mazes will always have exactly one entrance on the top row 
// -- you will know that you have exited the maze when you reach 
// the last row
// -- storage defined above is sufficient for static array (works 
// for maze1.txt, maze2.txt, maze3.txt) but we still check for 
// sufficient storage
// -- can only move one unit up, down, left, or right per move
//
// Implementation details:
// 
// The code:
// -- refers to directions north (up), south (down), west (left),
// and east (right)
// -- stores the maze data from the user-inputted maze file
// in a static 2D array 'a', where a[i][j] = 1 if (i,j) is
// the coordinate of a wall in the maze, and a[i][j] = 0 
// otherwise
// -- stores the coordinates of where to step in the solution via 
// the right hand rule in the 2D integer array 'soln', which is later
// used to write to the user-inputted solution file
// -- finds the (unique) maze entrance in row zero, storing this
// coordinate as the first position in 'soln'
// -- then moves one unit south, because this is the only possible next
// move without hitting a wall (otherwise there would be more than one
// maze entrance)
// -- records 'S' (for south) as the direction just moved
// -- then has different priorities of which direction to move in,
// if possible (with walls/maze boundary), depending on the direction 
// just moved. the following logic is used to obey the right hand rule:
// ---- if you just moved south ('S'), your priorities are: 'W' (west),
// 'S' (south), 'E' (east), 'N' (north) for which direction to go next.
// That is, go one unit west if you can, try going one unit south if you
// cannot go west, then try going one unit east if you can't go west or
// south, and finally go north as a last resort
// ---- if you just moved west ('W'), your priorities are: 'N', 'W',
// 'S','E'
// ---- if you just moved north ('N'), your priorities are 'E', 'N',
// 'W','S'
// ---- if you just moved east ('E'), your priorities are: 'S','E',
// 'N','W'
//
//
// BEGIN HELPER FUNCTIONS ////

tuple <int**,int> update_soln(tuple <int**,int> soln_tuple,
    int current_row, int current_col){
    // Helper function 'update_soln' records the move just
    // made by updating the solution data (contained in 
    // 'soln_tuple')
    int** modified = get<0>(soln_tuple);
    int soln_index = get<1>(soln_tuple);
    soln_index += 1;
    modified[soln_index][0] =current_row;
    modified[soln_index][1] =current_col;
    auto t = make_tuple(modified,soln_index);
    return t;
}

bool try_going_north(int current_row, int current_col,
    int a[num_rows_avail][num_cols_avail]){
    // Boolean function 'try_going_north' tests whether
    // it is possible to move one unit north given your
    // current position
    bool can_go_north = current_row - 1 >= 0 and
    a[current_row-1][current_col] == 0;
    return can_go_north;
}

bool try_going_south(int current_row, int current_col,
int a[num_rows_avail][num_cols_avail], int num_rows){
    // Boolean function 'try_going_north' tests whether
    // it is possible to move one unit south given your
    // current position
    bool can_go_south = current_row + 1 <= num_rows and 
    a[current_row+1][current_col] == 0;
    return can_go_south;
}

bool try_going_east(int current_row, int current_col,
    int a[num_rows_avail][num_cols_avail], int num_cols){
    // Similarly, boolean function 'try_going_east' checks whether
    // it is possible to move one unit east without hitting
    // a wall/maze boundary
    bool can_go_east = current_col + 1 <= num_cols and
            a[current_row][current_col+1] == 0;
    return can_go_east;
}

bool try_going_west(int current_row, int current_col,
int a[num_rows_avail][num_cols_avail]){
    // Boolean function 'try_going_west' checks whether
    // it is possible to move one unit west without hitting
    // a wall/maze boundary
    bool can_go_west = current_col - 1 >= 0 and 
    a[current_row][current_col-1] == 0;
    return can_go_west;
}
// END HELPER FUNCTIONS ////

int main(int argc, char *argv[]){
    // Initializes the number of rows/columns in the maze to the
    // maximum possible given storage constraints
    int num_rows = num_rows_avail;
    int num_cols = num_cols_avail;
    // Provides Usage statement for user
    if (argc != 3){
        cout << "Usage: "<< endl;
        cout << " " << argv[0] << " <maze file> <solution file>" << endl;
        return 0;
    }
    string mazefile = string(argv[1]);
    string solnfile = string(argv[2]);
    // Reading the maze file
    ifstream f(mazefile);
    // Writing to the solution file
    ofstream g(solnfile);

    // Creates static 2D data array 'a' for storing maze data
    int a[num_rows_avail][num_cols_avail];
    // Creates the 2D data array for storing solution data
    int **soln;
    // Each row in 'soln' is a two-integer array containing the
    // row and column of a move/position in the solution
    soln = new int* [num_rows_avail*num_cols_avail];
    //  Initializes soln array to all -1s
    for (int i = 0; i < num_rows_avail*num_cols_avail; i++){
        soln[i] = new int[2];
        soln[i][0] = -1;
        soln[i][1] = -1;
    }
    // Initialzes maze data array to all -1s (instructions were to
    // initialize all of the values in the array to a constant value, 
    // later changing each value corresponding to a maze location
    // to indicate the presence of a wall
    for(int i=0; i <num_rows_avail;i++){
        for(int j=0;j<num_cols_avail;j++){
            a[i][j] = -1;
        }
    }
    // Tries opening the maze file
    if (f.is_open()){
        int r;
        int c;
        f >> r >> c;
        // First line contains the number of rows/cols in the maze
        num_rows = r;
        num_cols= c;
        if(num_rows>num_rows_avail or num_cols>num_cols_avail){
            cout << "Not enough storage available" << endl;
            // Quits the program
            return 0;
        }
        // Initializes part of maze data array we actually
        // use (i.e. positions in the maze) to all zeros
        for (int i = 0; i < num_rows; i++){
            for (int j = 0; j < num_cols; j++){
                a[i][j] = 0;
            }
        }
        int row;
        int col;
        // Reads the row/cols in the maze data file corresponding
        // to wall locations
        while(f >> row >> col){
            // Again checks if a row/col of a wall corresponds to
            // a location beyond the storage space allocated
            if (row > num_rows or col > num_cols){
                cout << "Not enough storage available" << endl;
                // Quit the program
                return 0;
            }
            // Sets the wall positions to 1
            a[row][col] = 1;
        }
        f.close();
    }        
    else{
        cout << "Failed to open file: " << argv[1] << endl;
        // Quits the program
        return 0;
    }
    int entrance = -1;
    // Finds the maze entrance
    for(int j = 0; j<num_rows_avail;j++){
        if(a[0][j] == 0){
            entrance = j;
            // First position in the solution is the
            // entrance in row 0
            soln[0][0] = 0;
            soln[0][1] = entrance;
            // Next position in solution must be one step down,
            // because we are assuming only one entrance
            // in the first row
            soln[1][0] = 1;
            soln[1][1] = entrance;
        }
    }
    // Gets position data after the first move is made
    int current_row = 1;
    int current_col = soln[1][1];
    int soln_index = 1;
    // Creates the solution tuple, whose first entry is 'soln' (which
    // itself is a tuple) and whose second entry is the index of the
    // last non-trivial entry recorded (i.e. excluding -1s)
    tuple <int**,int> soln_tuple = make_tuple(soln,soln_index);
    // Keeps track of the last direction moved, which is south at
    // this point
    char direction = 'S';
    bool found_exit = false;
    
    while(found_exit == false){
        // Checks if you have exited the maze, breaking out of the
        // 'while' loop if so
        if(current_row == num_rows - 1){
            found_exit = true;
            break;
        }
        // Determines where to move next and which directions to try first
        // depending on the direction just moved
        switch(direction){
            // If you just moved south, try going west (if possible),
            // then south (if you cant go west), then east (if you 
            // cannot go west/south), and only move north as a last
            // resort
        case 'S':{
            // Attempts to move west
            bool can_go_west = try_going_west(current_row,current_col,a);
            if(can_go_west){
                current_col -= 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'W';
                break;
            }
            // Attempts to move south
            bool can_go_south = try_going_south(current_row,current_col,a,
                num_rows);
            if(can_go_south){
                current_row += 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'S';
                break;
            }
            // Attempts to move east
            bool can_go_east = try_going_east(current_row,current_col,a,
                num_cols);
            if(can_go_east){
                current_col += 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'E';
                break;
            }
            // Forced to go north. This is an option, for otherwise
            // you would be surrounded by 4 walls in the maze.
            current_row -= 1;
            soln_tuple = update_soln(soln_tuple,current_row,current_col);
            direction = 'N';
            break;
        }
        // If you just moved west:
        case 'W':{
            // First attempts to go north
            bool can_go_north = try_going_north(current_row,current_col,a);
            if(can_go_north){
                current_row -= 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'N';
                break;
            }
            // Then attempts to go west
            bool can_go_west = try_going_west(current_row,current_col,a);
            if(can_go_west){
                current_col -= 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'W';
                break;
            }
            // Then attempts to go south
            bool can_go_south = try_going_south(current_row,current_col,a,
            num_rows);
            if(can_go_south){
                current_row += 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'S';
                break;
            }
            // Forced to go east
            current_col += 1;
            soln_tuple = update_soln(soln_tuple,current_row,current_col);
            direction ='E';
            break;
        }
        // If you just moved north:
        case 'N':{
            // First tries going east
            bool can_go_east = try_going_east(current_row,current_col,a,
            num_cols);
            if(can_go_east){
                current_col += 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'E';
                break;
            }
            // Then tries going north
            bool can_go_north = try_going_north(current_row,current_col,a);
            if(can_go_north){
                current_row -= 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'N';
                break;
            }
            // Then tries going west
            bool can_go_west = try_going_west(current_row,current_col,a);
            if(can_go_west){
                current_col -= 1;
                soln_tuple = update_soln(soln_tuple,current_row,current_col);
                direction = 'W';
                break;
            }
            // Otherwise goes south
            current_row += 1;
            soln_tuple = update_soln(soln_tuple,current_row,current_col);
            direction = 'S';
            break;
        }
        // If you just moved east:
        default:
        // First tries going south
        bool can_go_south = try_going_south(current_row,current_col,a,
        num_rows);
        if(can_go_south){
            current_row += 1;
            soln_tuple = update_soln(soln_tuple,current_row,current_col);
            direction = 'S';
            break;
        }
        // Then tries going east
        bool can_go_east = try_going_east(current_row,current_col,a,
        num_cols);
        if(can_go_east){
            current_col += 1;
            soln_tuple = update_soln(soln_tuple,current_row,current_col);
            direction = 'E';
            break;
        }
        // Then attempts moving north
        bool can_go_north = try_going_north(current_row,current_col,a);
        if(can_go_north){
            current_row -= 1;
            soln_tuple = update_soln(soln_tuple,current_row,current_col);
            direction = 'N';
            break;
        }
        // Otherwise goes west
        current_col -= 1;
        soln_tuple = update_soln(soln_tuple,current_row,current_col);
        direction = 'W';
        }
    }
    // The 'while' loop was exited, so a solution has been found
    soln = get<0>(soln_tuple);
    // Writes the answer to the user-inputted solution file (referred
    // to previously as 'g')
    if (g.is_open()){
        for(int i = 0; i < num_rows*num_cols; i++){    
            // Finds the non-trivial entries of the solution array
            if(soln[i][0] != -1){
                // Writes those entries to the file
                g << soln[i][0] << " " << soln[i][1] << endl;
            }
        }
        g.close();
    }
    else{
        cout << "Failed to create file: " << argv[2] << endl;
        // Quits the program
        return 0;
    }
    // Finally exits the program
    return 0;
}
