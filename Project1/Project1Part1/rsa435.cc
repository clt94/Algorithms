// Connor Taylor
// Algorithms Project 1 Part 1: Using bigInt435 written by Matt McCuchen, implement 
// Fermat Test, generate two prime integers p and q greater than or equal to 512 bits, 
// and save them to p_q.txt. Then use the extended Euclidean algorithm to generate two
// keys (e, n) and (d, n) where n = p * q, and have e and n saved to e_n.txt and 
// d and n savec to d_n.txt.

// Standard libraries
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>

// 'BigIntegerLibrary.hh' includes all of the library headers.
#include "BigIntegerLibrary.hh"

int openFile(const std::string&, std::ofstream&);

BigUnsigned getPrime(const int);

BigUnsigned getRandInt(const int);

BigUnsigned getCoPrime(const BigUnsigned&);

// Size of numbers generated by random
const int SIZE = 200;

int main()
{
	/* The library throws 'const char *' error messages when things go
	 * wrong.  It's a good idea to catch them using a 'try' block like this
	 * one.  Your C++ compiler might need a command-line option to compile
	 * code that uses exceptions. */
	try
	{
		std::ofstream file;
		srand(time(0));
		
		// number of trials for Fermat's Little Theorem
		int k = 2;
		
		BigUnsigned p = getPrime(k);
		std::cout << "p: " << p << "\n";
		
		BigUnsigned q = getPrime(k);		
		std::cout << "q: " << q << "\n";
		
		// opens p_q.txt
		if(openFile("p_q.txt", file) == -1)
			return -1;
		// write p and q to p_q.txt then close
		file << p << '\n' << q << '\n';
		file.close();
		
		BigUnsigned n = p * q;
		std::cout << "n: " << n << "\n";
		
		BigUnsigned phi = (p - BigUnsigned(1)) * (q - BigUnsigned(1));
		std::cout << "phi: " << phi << '\n';
		
		BigUnsigned e = getCoPrime(phi);
		std::cout << "e: " << e << "\n";
		
		// opens e_n.txt
		if(openFile("e_n.txt", file) == -1)
			return -1;
		// writes e and n to e_n.txt then close
		file << e << '\n' << n << '\n';
		file.close();
		
		BigUnsigned d = modinv(e, phi);
		std::cout << "d: " << d << "\n";
		
		// opens d_n.txt
		if(openFile("d_n.txt", file) == -1)
			return -1;
		// writes d and n to d_n.txt then close
		file << d << '\n' << n << '\n';
		file.close();
		
	} catch(char const* err) 
	{
		std::cout << "The library threw an exception:\n" << err << std::endl;
	}

	return 0;
}

int openFile(const std::string& fileName, std::ofstream& file)
{
	// opens file in append mode
	file.open(fileName);
	// checks if file is open, if not exits program
	if (!file.is_open())
	{
		std::cout << fileName << " failed to open. Ending program." << std::endl;
		return -1;
	}
	
	return 0;
}

BigUnsigned getPrime(const int k)
{
	BigUnsigned p;
	BigUnsigned a;
	int count = 0;
	while(true)
	{
		// generates a random number, n, to check as prime, if count is 0
		if(count == 0)
		{
			p = BigUnsigned(1);
			// excludes any even, numbers divisible by 5, or 1
			while(p % 2 == 0 || p % 5 == 0 || p == 1)
				p = getRandInt(SIZE);
		}
			
		// generates a random number, a, less than p and greater then 2
		a = BigUnsigned(1);
		while(a > p || a < 2)
			a = getRandInt(rand() % SIZE);
			
		// check if prime
		if(modexp(a, p - 1, p) == 1)
		{
			// passed specified number of trials
			if(count == k)
				return p;
			// iterates number of tests if passed
			else 
				++count;
		}
		// resets number of tests if failed
		else 
			count = 0;
	}
}

BigUnsigned getRandInt(const int size)
{
	BigUnsigned num = BigUnsigned(1);
	// iterates to generate large number according to size
	for (int i = 0; i < size; ++i) 
		num = num * 10 + rand();
	
	return num;
}

BigUnsigned getCoPrime(const BigUnsigned& phi)
{
	BigUnsigned e;
	while(true)
	{
		// generates random number if not coprime
		e = getRandInt(rand() % (SIZE * 2));
		if(gcd(e, phi) == 1)
			return e;
	}
}
