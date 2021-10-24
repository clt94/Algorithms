#include <iostream>
#include <fstream>
#include <sstream>

void getEnergy(int**, int**, const int, const int);

void verEnergy(int**, int**, const int, const int);

void horEnergy(int**, int**, const int, const int);

void verSeam(int**, int**, const int, const int);

void horSeam(int**, int**, const int, const int);

void writeFile(int**, const int, const int, std::string&);

void deleteArray(int**, const int);

int main(int argc, char** argv) 
{
    if (argc != 4) 
	{
        std::cout << "Invalid number of arguments. Format should ./a <file file name> <vertical seams> <horizontal seams>";
    }
    int ver = atoi(argv[2]);
    int hor = atoi(argv[3]);
    std::string fileName = argv[1];
    std::string outFileName = fileName.substr(0, fileName.size() - 4);

    int numRows = 0; 
	int numCols = 0;
	int greyScaleMax = 255;
    std::ifstream file(fileName.c_str());
    std::stringstream ss;
    std::string input = "";
    getline(file, input);
    getline(file, input);
    ss << file.rdbuf();
    ss >> numCols >> numRows;
    ss >> greyScaleMax;

    // GreyScale Array used to store the greyscale values of the original image
    int** pgm = new int*[numRows];
    for (int i = 0; i < numRows; ++i) 
	{
        pgm[i] = new int[numCols];
		for (int j = 0; j < numCols; ++j) 
		{
            ss >> pgm[i][j];
        }
    }

    // close input file for efficiency
    file.close();

    //Carry our the vertical seam operation while ver is greater than zero
    for(int i = 0; i < ver; ++i) 
	{
        //create new array to store energy in
        int** energyPgm = new int*[numRows];
        int** cumEnergyPgm = new int*[numRows];
        for (int j = 0; j < numRows; ++j) 
		{
            energyPgm[j] = new int[numCols];
			cumEnergyPgm[j] = new int[numCols];
        }

        //calculate the energy array
        getEnergy(pgm, energyPgm, numRows, numCols);
        //calculate the cumulative vertical energies
        verEnergy(energyPgm, cumEnergyPgm, numRows, numCols);
        //flag positions of values to be deleted
        verSeam(pgm, cumEnergyPgm, numRows, numCols);
        //decrement the numColsumn number
        --numCols;

        //delete memory associated with cumEnergyPgm and energyPgm
		deleteArray(energyPgm, numRows);
		deleteArray(cumEnergyPgm, numRows);
    }

    //Carry out the horizontal seam operations while hor is greater than zero
    for(int i = 0; i < hor; ++i) 
	{
        //create new array to store energy in
        int** energyPgm = new int*[numRows];
		int** cumEnergyPgm = new int*[numRows];
        for (int j = 0; j < numRows; ++j) 
		{
            energyPgm[j] = new int[numCols];
			cumEnergyPgm[j] = new int[numCols];
        }

        //calculate the energy array
        getEnergy(pgm, energyPgm, numRows, numCols);
        //calculate the cumulative horizontal energies
        horEnergy(energyPgm, cumEnergyPgm, numRows, numCols);
        //flag positions of values to be deleted
        horSeam(pgm, cumEnergyPgm, numRows, numCols);
        //decrement the row number
        --numRows;

        //delete memory associated with cumEnergyPgm and energyPgm
		deleteArray(energyPgm, numRows);
		deleteArray(cumEnergyPgm, numRows);
    }
	
	

    //write pgm to a file
    writeFile(pgm, numRows, numCols, outFileName);

    //delete memory associated with pgm
    for (int i = 0; i < numRows; ++i) 
	{
        delete[] pgm[i];
    }
    delete[] pgm;

    return 0;
}

void getEnergy(int** pgm, int** energyPgm, const int numRows, const int numCols) 
{
    //loops through numRows
    for(int i = 0; i < numRows; ++i) 
	{
        for(int j = 0; j < numCols; ++j) 
		{
            //up and left
            if(i == 0 && j == 0) 
			{
                //only do right and down
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i + 1][j])) + abs((pgm[i][j]) - (pgm[i][j + 1])));
            }
            //up and right
            else if((i == 0) && (j == (numCols - 1)))
			{
                //only do left and down
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i][j - 1])) + abs((pgm[i][j]) - (pgm[i + 1][j])));
            }
            //down and left
            else if((i == (numRows - 1)) && (j == 0)) 
			{
                //only do up and right
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i - 1][j])) + abs((pgm[i][j]) - (pgm[i][j + 1])));
            }
            //down and right
            else if((i == (numRows - 1)) && (j == (numCols - 1)))
			{
                //only do up and left
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i - 1][j])) + abs((pgm[i][j]) - (pgm[i][j - 1])));
            }
            //for top row
            else if(i == 0) 
			{
                //right, left and bottom
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i][j + 1])) + abs((pgm[i][j]) - (pgm[i][j - 1])) +
                                   abs((pgm[i][j]) - (pgm[i + 1][j])));
            }
            //for bottom row
            else if(i == (numRows - 1))
			{
                //top, left and right
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i][j + 1])) + abs((pgm[i][j]) - (pgm[i][j - 1])) +
                                   abs((pgm[i][j]) - (pgm[i - 1][j])));
            }
            //for first numColsumn
            else if(j == 0) 
			{
                //right, top and bottom
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i][j + 1])) + abs((pgm[i][j]) - (pgm[i - 1][j])) +
                                   abs((pgm[i][j]) - (pgm[i + 1][j])));
            }
            else if(j == (numCols - 1))
			{
                //left, top and bottom
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i][j - 1])) + abs((pgm[i][j]) - (pgm[i + 1][j])) +
                                   abs((pgm[i][j]) - (pgm[i - 1][j])));
            } 
			else 
			{
                energyPgm[i][j] = (abs((pgm[i][j]) - (pgm[i - 1][j])) + abs((pgm[i][j]) - (pgm[i + 1][j])) +
                                   abs((pgm[i][j]) - (pgm[i][j - 1])) + abs((pgm[i][j]) - (pgm[i][j + 1])));
            }
        }
    }
}

void verEnergy(int** energyPgm, int** cumEnergyPgm, const int numRows, const int numCols) 
{
    //loop through the numRows
    for(int i = 0; i < numRows; ++i) 
	{
        //loop through the numCols
        for(int j = 0; j < numCols; ++j)
		{
            //do this only for the first row
            if(i == 0)
			{
                //copy over values
                cumEnergyPgm[i][j] = energyPgm[i][j];
            } 
			else 
			{
                //first numColsumn
                if(j == 0) 
				{
                    //ignore out of bounds
                    cumEnergyPgm[i][j] = energyPgm[i][j] + std::min(cumEnergyPgm[i - 1][j], cumEnergyPgm[i - 1][j + 1]);
                } 
				else if(j == (numCols - 1)) 
				{
                    //ignore out of bounds
                    cumEnergyPgm[i][j] = energyPgm[i][j] + std::min(cumEnergyPgm[i - 1][j - 1], cumEnergyPgm[i - 1][j]);
                } 
				else 
				{
                    //find the min energy to add to the running total
                    cumEnergyPgm[i][j] = energyPgm[i][j] + 
										  std::min(std::min(cumEnergyPgm[i - 1][j - 1], cumEnergyPgm[i - 1][j]), cumEnergyPgm[i - 1][j + 1]);
                }
            }
        }
    }
}

void horEnergy(int** energyPgm, int** cumEnergyPgm, const int numRows, const int numCols) 
{
    //loop through the numCols
    for(int j = 0; j < numCols; ++j) 
	{
        //loop through the numRows
        for(int i = 0; i < numRows; ++i) 
		{
            //do this only for the first numColsumn
            if(j==0)
			{
                //copy over values
                cumEnergyPgm[i][j] = energyPgm[i][j];
            } 
			else 
			{
                //top row
                if(i == 0) 
				{
                    //ignore out of bounds
                    cumEnergyPgm[i][j] = energyPgm[i][j] + std::min(cumEnergyPgm[i][j - 1], cumEnergyPgm[i + 1][j - 1]);
                }
				else if(i==(numRows-1)) 
				{
                    //ignore out of bounds
                    cumEnergyPgm[i][j] = energyPgm[i][j] + std::min(cumEnergyPgm[i][j - 1], cumEnergyPgm[i - 1][j - 1]);
                } 
				else
				{
                    //find the min energy to add to the running total
                    cumEnergyPgm[i][j] = energyPgm[i][j] + std::min(std::min(cumEnergyPgm[i - 1][j - 1], cumEnergyPgm[i][j - 1]), cumEnergyPgm[i + 1][j - 1]);
                }
            }
        }
    }
}

void verSeam(int** pgm, int** cumEnergyPgm, const int numRows, const int numCols) 
{
    int minIndex = 0;
    //set minimum to the first element in the last row
    int min = cumEnergyPgm[numRows - 1][minIndex];

    //loop through last row and find index of the minimum
    for(int i = 0; i < numCols; ++i) 
	{
        if(cumEnergyPgm[numRows - 1][i] <= min) 
		{
            minIndex = i;
            min = cumEnergyPgm[numRows - 1][i];
        }
    }

    //delete element and shift
    for(int i = minIndex; i < (numCols - 1); ++i)
	{
        //copy the left element of the last row
        pgm[numRows - 1][i] = pgm[numRows - 1][i + 1];
    }

    //loop through remaining numRows
    for(int i = (numRows - 2); i >= 0; --i) 
	{
        //base case for edge numRows
        if(minIndex == 0)
		{
            //check the minium of directly above or to the right
            if(cumEnergyPgm[i][minIndex] < cumEnergyPgm[i][minIndex + 1])
			{
                //return above
                minIndex = minIndex;
            } 
			else if(cumEnergyPgm[i][minIndex + 1] < cumEnergyPgm[i][minIndex])
			{
                //return right
                minIndex = minIndex + 1;
            }
        }
        // right edge
        else if(minIndex == (numCols - 1))
		{
            //only check the minium of directly above or to the right
            if(cumEnergyPgm[i][minIndex] < cumEnergyPgm[i][minIndex - 1])
			{
                //return above
                minIndex = minIndex;
            }
			else if(cumEnergyPgm[i][minIndex - 1] < cumEnergyPgm[i][minIndex])
			{
                //return left
                minIndex = minIndex - 1;
            }
        }
        //not first or last numColsumn case
        // find min of above left and right
        else 
		{
            // check middle as min
            if((cumEnergyPgm[i][minIndex] < cumEnergyPgm[i][minIndex - 1]) && (cumEnergyPgm[i][minIndex] < cumEnergyPgm[i][minIndex + 1]))
			{
                //return middle
                minIndex = minIndex;
            } // check the right as min
            else if((cumEnergyPgm[i][minIndex + 1] < cumEnergyPgm[i][minIndex]) && (cumEnergyPgm[i][minIndex + 1] < cumEnergyPgm[i][minIndex - 1]))
			{
                //return right
                minIndex = minIndex + 1;
            } // check left as min
            else if((cumEnergyPgm[i][minIndex - 1] < cumEnergyPgm[i][minIndex]) && (cumEnergyPgm[i][minIndex - 1] < cumEnergyPgm[i][minIndex + 1]))
			{
                //return left
                minIndex = minIndex - 1;
            } // if both left and right are equal and less than middle
            else if((cumEnergyPgm[i][minIndex - 1] < cumEnergyPgm[i][minIndex]) && (cumEnergyPgm[i][minIndex + 1] < cumEnergyPgm[i][minIndex]) && 
				    (cumEnergyPgm[i][minIndex + 1] == cumEnergyPgm[i][minIndex - 1]))
			{
                // by default return the left
                minIndex = minIndex - 1;
            }
        }

        // delete minIndex value and shift all the values over
        for(int j = minIndex; j < (numCols - 1); ++j)
		{
            //copy over the left element
            pgm[i][j] = pgm[i][j + 1];
        }
    }
}

void horSeam(int** pgm, int** cumEnergyPgm, const int numRows, const int numCols) 
{
    int minIndex = numCols - 1;
    int min = 999999;

    //loop through last numCols and find index of the minimum
    for(int i = 0; i < numRows; ++i) 
	{
        if(cumEnergyPgm[i][numCols - 1] < min) 
		{
            minIndex = i;
            min = cumEnergyPgm[i][numCols - 1];
        }
    }

    for(int i = minIndex; i < (numRows - 1); ++i)
	{
        //copy the left element of the last row
        pgm[i][numCols - 1] = pgm[i + 1][numCols - 1];
    }

    //loop through the remaining numCols
    for(int i = (numCols - 2); i >= 0; --i)
	{
        //base case for first or last row
        //first row
        if(minIndex==0)
		{
            //only check the minium of directly left and to the bottom
            //check left as min
            if(cumEnergyPgm[minIndex][i] < cumEnergyPgm[minIndex + 1][i])
			{
                //return left
                minIndex = minIndex;
            } 
			else if(cumEnergyPgm[minIndex + 1][i] < cumEnergyPgm[minIndex][i])
			{
                //return bottom
                minIndex = minIndex + 1;
            }
        }
        // right
        else if(minIndex == (numRows - 1))
		{
            //check left as min
            if(cumEnergyPgm[minIndex][i] < cumEnergyPgm[minIndex - 1][i])
			{
                //return left
                minIndex = minIndex;
            } 
			else if(cumEnergyPgm[minIndex - 1][i] < cumEnergyPgm[minIndex][i])
			{
                //return top
                minIndex = minIndex - 1;
            }
        }
        //not first or last row
        // find min of top left and bottom
        else 
		{
            // check middle as min
            if((cumEnergyPgm[minIndex][i] < cumEnergyPgm[minIndex - 1][i]) && (cumEnergyPgm[minIndex][i] < cumEnergyPgm[minIndex + 1][i]))
			{
                //return left
                minIndex = minIndex;
            } // check the top as min
            else if((cumEnergyPgm[minIndex - 1][i] < cumEnergyPgm[minIndex][i]) && (cumEnergyPgm[minIndex - 1][i] < cumEnergyPgm[minIndex + 1][i]))
			{
                //return top
                minIndex = minIndex - 1;
            } // check bottom as min
            else if((cumEnergyPgm[minIndex + 1][i] < cumEnergyPgm[minIndex][i]) && (cumEnergyPgm[minIndex + 1][i] < cumEnergyPgm[minIndex - 1][i]))
			{
                //return bottom
                minIndex = minIndex + 1;
            }
            // if both top and bottom are equal and less than left
            else if((cumEnergyPgm[minIndex + 1][i] < cumEnergyPgm[minIndex][i]) && (cumEnergyPgm[minIndex - 1][i] < cumEnergyPgm[minIndex][i]) && 
					(cumEnergyPgm[minIndex + 1][i] == cumEnergyPgm[minIndex - 1][i]))
			{
                // by default return the top
                minIndex = minIndex - 1;
            }
        }
        // delete minIndex value and shift
        for(int j = minIndex; j < (numRows - 1); ++j)
		{
            //copy over the left element
            pgm[j][i] = pgm[j + 1][i];
        }
    }
}

//Writes the seam carved image to a .file file for viewing
void writeFile(int** array, const int numRows, const int numCols, std::string &outFileName) 
{
    std::ofstream outFile;
    outFileName.append("_processed.pgm");
    outFile.open(outFileName.c_str());
    if (!outFile.is_open()) 
	{
        std::cout << "Can't open output file"  << outFileName << std::endl;
        return;
    }

    // write the header of the file file
    outFile << "P2\n" << "#Output file for " << outFileName << "\n" << numCols << " " << numRows << "\n255\n";

    //writes the contents of the array to a .file file for viewing
    for(int i = 0; i < numRows; ++i) 
	{
        //loops through numRows and numCols
        for(int j = 0; j < numCols - 1; ++j) 
		{
            outFile << array[i][j] << " ";
        }
        //last numColsumn case appends new line
        outFile << array[i][numCols - 1] << std::endl;
    }
    outFile.close();
 }
 
void deleteArray(int** array, int numRows)
{
	for (int i = 0; i < numRows; ++i) 
	{
		delete[] array[i];
	}
	
	delete[] array;
}