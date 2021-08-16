/*
 * InteractionArc.cpp
 *
 *  Created on: Apr 20, 2014
 *      Author: psz
 */

#include "../include/InteractionArc.h"

InteractionArc::InteractionArc() {
	init();
}

InteractionArc::InteractionArc(int _start, int _end, int _score, int _factor) {
	init();
	start = _start;
	end = _end;
	score = _score;
	eff_score = score;
	factor = _factor;
}

void InteractionArc::init() {
	start = -1;
	end = -1;
	genomic_start = -1;
	genomic_end = -1;
	score = 0;
	eff_score = 0;
	factor = -1;
}

void InteractionArc::toFile(std::ofstream & file) {
	std::stringstream out;
	this->toStringStream(out);
	file << out.str();
}

void InteractionArc::fromFile(std::ifstream & file) {
	std::string full_input(std::istreambuf_iterator<char>{file}, {});
	std::stringstream in_stream(full_input);

	this->fromStringStream(in_stream);
}

void InteractionArc::toStringStream(std::stringstream & out) {
	out << start << ' ' << end << ' ' << genomic_start << ' ' << genomic_end << ' ';
	out << score << ' ' << eff_score << ' ' << factor << '\n';
}

void InteractionArc::fromStringStream(std::stringstream & in) {
	in >> start >> end >> genomic_start >> genomic_end >> score >> eff_score >> factor;
}

int InteractionArc::length() {
	return end - start + 1;
}

void InteractionArc::print() {
	printf("%d %d %d (eff_sc: %d, factor: %d)\n", start, end, score, eff_score, factor);
}
