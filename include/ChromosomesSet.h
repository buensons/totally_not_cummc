/*
 * ChromosomesSet.h
 *
 *  Created on: May 25, 2014
 *      Author: psz
 */

#ifndef CHROMOSOMESSET_H_
#define CHROMOSOMESSET_H_

#include <string.h>
#include <vector>
#include <fstream> 
#include <sstream>
#include <iostream>
#include "Chromosome.h"

class ChromosomesSet {
public:
	ChromosomesSet();

	void print();

	void add(std::map<std::string, Chromosome> chr);
	void add(std::map<std::string, Chromosome> chr, string desc);

	void toFile(string filename);
	void fromFile(string filename);
	void fromStringStream(std::stringstream & in);

	std::vector<std::map<std::string, Chromosome> > chromosome;
	std::vector<string> desc;

};

#endif /* CHROMOSOMESSET_H_ */
