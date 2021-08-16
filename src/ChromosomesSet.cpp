/*
 * ChromosomesSet.cpp
 *
 *  Created on: May 25, 2014
 *      Author: psz
 */

#include "../include/ChromosomesSet.h"

ChromosomesSet::ChromosomesSet() {
	// TODO Auto-generated constructor stub

}


void ChromosomesSet::print() {
	printf("Set size = %lu\n", chromosome.size());

	for (size_t i = 0; i < chromosome.size(); ++i) {
	//	printf("%d %s\n", chromosome[i].points.size(), desc[i].c_str());
	}
}

void ChromosomesSet::add(std::map<std::string, Chromosome> chr) {
	add(chr, "<no_desc>");
}

void ChromosomesSet::add(std::map<std::string, Chromosome> chr, string desc) {
	chromosome.push_back(chr);
	if (desc == "") desc = "none";
	this->desc.push_back(desc);
}


void ChromosomesSet::toFile(string filename) {
	std::stringstream out;
	std::ofstream file;
	file.exceptions(std::ofstream::badbit | std::ofstream::failbit);

	try {
		file.open(filename);
	} catch(const std::ofstream::failure &ex) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	out << chromosome.size() << '\n';

	for (size_t i=0; i < chromosome.size(); i++) {
		out << chromosome[i].size() << ' ' << desc[i] << '\n';

		for (auto el: chromosome[i]) {
			out << el.first << ' ' << el.second.points.size() << '\n';
			el.second.toStringStream(out);
		}
	}

	file << out.str();
	file.close();
}

void ChromosomesSet::fromFile(string filename) {
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

	fromStringStream(in_stream);
	
	file.close();
}

void ChromosomesSet::fromStringStream(std::stringstream & in) {
	chromosome.clear();

	char de[100], str_chr[10];
	int n, n_chr;
	in >> n;

	int pts;
	for (int i=0; i<n; i++) {
		in >> n_chr >> de;
		desc.push_back(std::string(de));

		std::map<std::string, Chromosome> map_chr;
		for (int j = 0; j < n_chr; ++j) {
			in >> str_chr >> pts;

			Chromosome chr;
			chr.fromStringStream(in, pts);

			map_chr[str_chr] = chr;
		}

		chromosome.push_back(map_chr);
	}
}
