/*
 * Cluster.h
 *
 *  Created on: Jun 1, 2014
 *      Author: psz
 */

#pragma once

#include <stdio.h>
#include <vector>
#include <fstream> 
#include <sstream>
#include <iostream>

#include "InteractionArc.h"
#include "../src/lib/common.h"

class Cluster {
public:
	Cluster();
	Cluster(int start, int end);

	void init();
	void print();

	void toFile(std::ofstream & file);
	void toFilePreviousFormat(std::ofstream & file);
	void fromFile(std::ifstream & file);

	void toStringStream(std::stringstream & out);
	void toStringStreamPreviousFileFormat(std::stringstream & out);
	void fromStringStream(std::stringstream & in);

	bool contains(int genomic_pos);	// check if a genomic position is contained in a given cluster

	vector3 pos;					// 3D position
	int genomic_pos;				// genomic position
	int start, end;					// genomic start and end position
	char orientation;

	int parent;
	int level;

	int base_start, base_end;

	std::vector<int> arcs;
	std::vector<int> siblings;

	std::vector<int> children;

	bool is_fixed;
	double dist_to_next;
};
