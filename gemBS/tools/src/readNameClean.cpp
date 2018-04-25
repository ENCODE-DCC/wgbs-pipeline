/*
 * main.cpp
 *
 *  Created on: 21 Out, 2015
 *      Author: marcos
 */
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

/**
 * \brief splits a string and leaves the results in a vector
 * \param s string to be split
 * \param delim character delimitate
 * \param elems vector of string elements to be fulfilled with the different fields
 */
vector<string> split(std::string s, char delim, vector<string> elems)
{
	std::stringstream ss(s);
	std::string item;

	while(getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	return elems;
}

/**
 * \brief Splits a string according to a given char seprator
 * \param s string to be split
 * \param delim character delimitate
 * \returs vector of elements
 */
vector<string> split(const string &s, char delim)
{
		vector<string> elems;
	    return split(s, delim, elems);
}

/**
 * \brief Remove all chars presents in charsToRemove from string s
 * \param str, string to treat
 * \param charsToRemove, vector of characters to remove
 */
void removeCharsFromString( string &str, char* charsToRemove )
{
   for ( unsigned int i = 0; i < strlen(charsToRemove); ++i )
   {
      str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
   }
}

/**
 * \brief main function
 * \brief Reads sam files and cleans from read name @ character
 */
int main(int argc, char *argv[])
{

	string line;

	//  Don't sync C++ and C I/O
	ios_base::sync_with_stdio(false);

	while(cin)
	{
		//1.2 Get line from standard input
	    getline(cin, line);

	    if(!line.empty())
	    {
	    	//1.3 Split in field and remove @
	    	vector<string> vFields = split(line, '\t');
	    	string firstField = vFields.at(0);

	    	if (firstField[0] != '@')
	    	{
	    		removeCharsFromString( firstField, "@" );
	    		vFields[0] = firstField;
	    	}

	    	//1.4 Print to standard output
	    	vector<string>::iterator final_iter = vFields.end();
	    	--final_iter;

	    	for(vector<string>::iterator it = vFields.begin(); it != vFields.end(); ++it)
	    	{
	    		printf("%s",it->c_str());

	    		if (it != final_iter)
	    		{
	    			printf("\t");
	    		}
	    	}
	    	printf("\n");
	    }
	}

	return 0;
}
