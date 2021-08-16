/*
 * BedRegions.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: psz
 */

#include "../include/BedRegions.h"

BedRegions::BedRegions() {
}

void BedRegions::fromFile(std::string filename) {
	std::ifstream file;
	file.exceptions(std::ifstream::badbit | std::ifstream::failbit);

	try {
		file.open(filename);
	} catch(const std::ifstream::failure &ex) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	std::string full_input(std::istreambuf_iterator<char>{file}, {});
	std::stringstream in_stream(full_input);

	std::string chr;
	int start, end;

	while (!in_stream.eof()) {
		in_stream >> chr >> start >> end;

		if (in_stream.fail()) {
            in_stream.clear();
            in_stream.ignore(100, '\n');
            continue;
        }

		in_stream.ignore(100, '\n');	// read to the end of line

		BedRegion b (chr, start, end);
		regions.push_back(b);
	}
	file.close();
}

void BedRegions::print() {
	printf("regions: %d\n", (int)regions.size());
	for (unsigned int i = 0; i < regions.size(); ++i) {
		printf("[%s] %d %d\n", regions[i].chr.c_str(), regions[i].start, regions[i].end);
	}
}

void BedRegions::addNewIntervals(std::string chr, int start, int end, int step) {

	if (step < 1) {
		printf("Step must be at least 1!\n");
		return ;
	}

	while (start < end) {
		BedRegion reg(chr, start, start);
		regions.push_back(reg);
		start += step;
	}
}
