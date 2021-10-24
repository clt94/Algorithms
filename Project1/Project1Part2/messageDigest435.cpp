/***
   prepared for CS435 Project 1 part 2
**/

#include <string.h>
#include <iostream>
#include <fstream>
#include "sha256.h"
#include "BigIntegerLibrary.hh"

 
int main(int argc, char *argv[])
{
	/*
	//demonstrating how sha256 works
	std::string input = "testing";
	std::string output1 = sha256(input);
	std::cout << "sha256('"<< input << "'):" << output1 << "\n";
   
	//demo bigInt works here
	BigUnsigned a = stringToBigUnsigned("124338907642982722929222542626327282");
	BigUnsigned b = stringToBigUnsigned("124338907642977775546469426263278643");
	std::cout << "big a = " <<a<<"\n";
	std::cout << "big b = " <<b<<"\n";
	std::cout << "big a*b = " <<a*b<<"\n";
	*/
	//Second part of your project starts here
	if (argc != 3 || (argv[1][0]!='s' && argv[1][0]!='v')) 
		std::cout << "wrong format! should be \"messageDirect435.exe s filename\"";
	else 
	{
		std::string filename = argv[2];
	  
		//read the file
		std::streampos begin,end;
		std::ifstream myfile (filename.c_str(), std::ios::binary);
		begin = myfile.tellg();
		myfile.seekg (0, std::ios::end);
		end = myfile.tellg();
		std::streampos size = end-begin;
		//std::cout << "size of the file: " << size << " bytes.\n"; //size of the file
	  
		myfile.seekg (0, std::ios::beg);
		char * memblock = new char[size];
		myfile.read (memblock, size); //read file; it's saved in the char array memblock
		myfile.close();
	  
		std::string copyOFfile = filename+".Copy"; 
		std::ofstream myfile2 (copyOFfile.c_str(), std::ios::binary);
		myfile2.write (memblock, size); //write to a file
		myfile2.close();
	  
		//std::cout << memblock;
		
		// takes file from commmand line, converts it to base 10
		BigUnsigned hash = BigUnsignedInABase(sha256(memblock), 16);
		std::string temp;
		BigUnsigned sig;
		// signs the document
		if (argv[1][0]=='s') 
		{
			// gets d from text document
			myfile.open("d_n.txt");
			if(!myfile.is_open())
				return -1;
			myfile >> temp;
			BigUnsigned d = stringToBigUnsigned(temp);
			temp = "";
			
			// gets n from text document
			myfile >> temp;
			BigUnsigned n = stringToBigUnsigned(temp);
			myfile.close();
			
			// signs the document and produces .signature file
			sig = modexp(hash, d, n);
			myfile2.open(filename+".signature");
			if(!myfile2.is_open())
				return -1;
			myfile2 << sig << '\n';
			myfile2.close();
		}
		// encrypts the document
		else 
		{
			// opens file and gets signed value
			myfile.open(filename);
			if(!myfile.is_open())
				return -1;
			myfile >> temp;
			myfile.close();
			sig = stringToBigUnsigned(temp);
			temp = "";
			
			// gets e from text file
			myfile.open("e_n.txt");
			if(!myfile.is_open())
				return -1;
			myfile >> temp;
			BigUnsigned e = stringToBigUnsigned(temp);
			temp = "";
			
			// gets n from text file
			myfile >> temp;
			BigUnsigned n = stringToBigUnsigned(temp);
			myfile.close();
			
			// gets encrpyted value from signature 
			BigUnsigned enc = modexp(sig, e, n);
			
			// checks that the document has been modified 
			if(enc != hash)
				std::cout << "Modified document\n";
			else
				std::cout << "Authentic document\n";
		}
		delete[] memblock;
	}
	return 0;
}