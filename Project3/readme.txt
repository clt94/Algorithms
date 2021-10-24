This program takes a .pgm file and removes a user specified amount of vertical and horizontal seams. The program then outputs a new .pgm file with the appropiate computations.

To compile:
g++ -o a main.cpp

To run:
./a <Your pgm file> <Number of vertical seams> <Number of horizontal seams>

Example:
./a Buchtel.pgm 50 20