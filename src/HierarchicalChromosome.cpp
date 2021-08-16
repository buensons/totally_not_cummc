/*
 * HierarchicalChromosomeMixed.cpp
 *
 *  Created on: Jun 2, 2014
 *      Author: psz
 */

#include "../include/HierarchicalChromosome.h"

HierarchicalChromosome::HierarchicalChromosome() {
	smooth_factor = -1;
}

void HierarchicalChromosome::print(int level) {
	printf("*** clusters: %d\n", (int)clusters.size());

	for (string chr: chrs) {
		printf(" [%s] root: %d, arcs: %d, curr_level size: %d\n", chr.c_str(), chr_root[chr], (int)arcs.arcs[chr].size(), (int)current_level[chr].size());
	}

	if (level >= 1) {
		printf("clusters\n");
		for (unsigned int i = 0; i < clusters.size(); ++i) {
			printf("[%u] ", i);
			clusters[i].print();
		}

		printf("current regions\n");
		for (string chr: chrs) {
			printf("%s\n", chr.c_str());
			for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
				printf("[%u] %d\n", i, current_level[chr][i]);
			}
		}

		//printf("arcs: (cnt=%d)\n", arcs.arcs.size());
		//arcs.print();
	}
}

vector3 HierarchicalChromosome::getCenter(bool current) {
	vector3 center(0.0f, 0.0f, 0.0f);
	if (current) {
		int cnt = 0;
		for (string chr: chrs) {
			for (unsigned int i=0; i<current_level[chr].size(); i++) center += clusters[current_level[chr][i]].pos;
			cnt += current_level[chr].size();
		}
		center /= cnt;
	}
	else {
		for (unsigned int i = 0; i < clusters.size(); ++i) center += clusters[i].pos;
		center /= clusters.size();
	}
	return center;
}

void HierarchicalChromosome::center(bool current) {
	vector3 center = getCenter(current);

	if (current) {
		for (string chr: chrs) {
			for (unsigned int i = 0; i < current_level[chr].size(); ++i) clusters[current_level[chr][i]].pos -= center;
		}
	}
	else {
		for (unsigned int i = 0; i < clusters.size(); ++i) clusters[i].pos -= center;
	}
}

map<string, float> HierarchicalChromosome::getAvgAnchorDistance() {
	map<string, float> map;
	float d, tot_dist = 0.0f;
	int tot_cnt;
	for (string chr: chrs) {
		d = 0.0f;
		for (unsigned int i = 1; i < current_level[chr].size(); ++i) {
			d += (clusters[current_level[chr][i]].pos - clusters[current_level[chr][i-1]].pos).length();
		}
		tot_dist += d;
		tot_cnt += current_level[chr].size();

		d /= current_level[chr].size();
		map[chr] = d;
	}
	map["genome"] = tot_dist / tot_cnt;
	return map;
}

void HierarchicalChromosome::scale(float factor) {
	vector3 center = getCenter();
	for (unsigned int i = 0; i < clusters.size(); ++i) {
		clusters[i].pos = center + (clusters[i].pos - center) * factor;
	}
	createCurrentLevelStructure();
}

//int HierarchicalChromosomeMixed::findClosestCluster(int genomic_position) {
//
//	levelDown();
//	levelDown();
//
//	int n = current_level.size();
//	if (n < 2) return n-1;	// -1 if no clusters, 0 if there is only one
//	if (genomic_position < clusters[current_level[0]].base_start || genomic_position > clusters[current_level[0]].base_start);
//
//
//	//printv(current_level)
//
//	int l = 0, r = n-1;
//	int mid = 0;
//	while (r>l+1) {
//		mid = (r+l) / 2;
//
//		printf("%d %d, mid=%d\n", l, r, mid);
//		if (genomic_position < mid) r = mid;
//		else l = mid;
//
//	}
//}

//void HierarchicalChromosomeMixed::align(const HierarchicalChromosomeMixed &hc) {
//	int i, j, n = current_level.size();
//	int steps = 20;
//	float ax, ay, az;
//	float score, err, best = 10e9;
//	matrix44 m, mbest;
//	map<int, int> map;	// map index from this->clusters to hc->clusters
//
//	center();
//
//
//	std::vector<vector3> pos;	// save original position
//	for (i = 0; i < n; ++i) pos.push_back(clusters[current_level[i]].pos);
//
//	for (i = 0; i < steps; ++i) {
//		ax = random(2.0f * 3.14f, false);
//		ay = random(2.0f * 3.14f, false);
//		az = random(2.0f * 3.14f, false);
//
//		// create rot matrix
//		m = RotateRadMatrix44('x', ax);
//		m = m * RotateRadMatrix44('y', ay);
//		m = m * RotateRadMatrix44('z', az);
//
//		// rotate
//		for (j = 0; j < n; ++j) clusters[current_level[i]].pos = clusters[current_level[i]].pos * m;
//
//		// calc score
//		score = 0.0f;
//		for (j = 0; j < n; ++j) {
//			//err = (clusters[current_level[i]].pos - hc.clusters[current_level[i]].pos).length();
//			//clusters[current_level[i]].pos = clusters[current_level[i]].pos * m;
//		}
//
//		// restore previous
//		for (j = 0; j < n; ++j) clusters[current_level[i]].pos = clusters[current_level[i]].pos * m;
//	}
//
//
//}

int HierarchicalChromosome::clusterIndexToSmoothedIndex(int ind) {
	if (smooth_factor == -1) return -1;
	return ind * (smooth_factor + 1);
}

void HierarchicalChromosome::printRegionsTree() {
	for (string chr: chrs) {
		printRegionsTree(chr);
	}
}

void HierarchicalChromosome::printRegionsTree(string chr) {
	printf("%s\n", chr.c_str());
	for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
		printRegionsTreeRecursive(current_level[chr][i]);
	}
}

void HierarchicalChromosome::printRegionsTreeRecursive(int region_index_curr, int margin) {

	for (int i = 0; i < margin; ++i) printf("   ");
	clusters[region_index_curr].print(); // print current

	// print subregions
	for (unsigned int i = 0; i < clusters[region_index_curr].children.size(); ++i) {
		printRegionsTreeRecursive(clusters[region_index_curr].children[i], margin+1);
	}
}


void HierarchicalChromosome::useTopLevel() {
	if (clusters.size() > 0) {
		current_level.clear();
		for (string chr: chrs) current_level[chr].push_back(chr_root[chr]);
		createCurrentLevelStructure();
	}
}

void HierarchicalChromosome::useLowestLevel() {
	for (int i = 0; i < 5; ++i) levelDown();
	createCurrentLevelStructure();
}

void HierarchicalChromosome::setLevel(int level) {
	useTopLevel();
	while (level--) levelDown();
}

void HierarchicalChromosome::levelDown() {
	std::vector<int> tmp;
	for (string chr: chrs) {

		for (unsigned int i = 0; i<current_level[chr].size(); ++i) {

			// if cluster has no children, then keep it
			if (clusters[current_level[chr][i]].children.size() == 0) {
				tmp.push_back(current_level[chr][i]);
			}
			else {
				for (unsigned int j = 0; j < clusters[current_level[chr][i]].children.size(); ++j) {
					tmp.push_back(clusters[current_level[chr][i]].children[j]);
				}
			}
		}

		current_level[chr] = tmp;
		tmp.clear();
	}
	createCurrentLevelStructure();
}


void HierarchicalChromosome::expandRegion(int start, int end, bool include_external) {
	//	std::vector<int> v;
	//
	//	for (int i = 0; i < current_level.size(); ++i) {
	//		std::vector<int> vc = expandRegion(current_level[i], start, end, include_external);
	//		if (vc.size() > 0) v.insert(v.end(), vc.begin(), vc.end());
	//		else v.insert(v.end(), current_level[i]);
	//	}
	//
	//	current_level = v;
	//	createCurrentLevelStructure();
}

std::vector<int> HierarchicalChromosome::expandRegion(int region_ind, int start, int end, bool include_external) {
	std::vector<int> v;
	//	int p;
	//	//printf("expand region %d\n", region_ind);
	//
	//	if (end < clusters[region_ind].start || start > clusters[region_ind].end) return v;  // not included
	//
	//	for (int i = 0; i < clusters[region_ind].children.size(); ++i) {
	//		p = clusters[region_ind].children[i];
	//		//printf("%d %d %d %d\n", end , clusters[p].start , start , clusters[p].end);
	//		if (end < clusters[p].start || start > clusters[p].end) {
	//			v.insert(v.end(), p);
	//			//printf("ngh\n");
	//			continue;
	//		}
	//
	//		if (clusters[p].children.size() == 0) {
	//			v.insert(v.end(), p);
	//			//printf("no child\n");
	//			continue;
	//		}
	//
	//		std::vector<int> vc = expandRegion(p, start, end, include_external);
	//		if (vc.size() > 0) v.insert(v.end(), vc.begin(), vc.end());
	//		else v.insert(v.end(), current_level[i]);
	//	}
	//
	//	/*
	//	for (int i = 0; i < clusters[region_ind].children.size(); ++i) {
	//		p = clusters[region_ind].children[i];
	//		printf("%d %d %d %d\n", end , clusters[p].start , start , clusters[p].end);
	//		if (end < clusters[p].start || start > clusters[p].end) {
	//			if (include_external) v.insert(v.end(), p);
	//			printf("ig\n");
	//			continue;
	//		}
	//
	//		if (clusters[p].children.size() == 0) {
	//			v.insert(v.end(), p);
	//			printf("no child\n");
	//			continue;
	//		}
	//
	//		printf("ok [%d] ", p);
	//
	//		std::vector<int> vc = expandRegion(p, start, end, include_external);
	//		printv(vc, true, true, "expanded children");
	//		v.insert(v.end(), vc.begin(), vc.end());
	//	}
	//	 */
	//	//printv(v, true, true, "v");
	return v;

	/*
	std::vector<int> tmp;
	for (int i = 0; i<current_level.size(); ++i) {
		if (clusters[current_level[i]].children.size() == 0) {
			tmp.push_back(current_level[i]);
		}
		else {
			for (int j = 0; j < clusters[current_level[i]].children.size(); ++j) {
				tmp.push_back(clusters[current_level[i]].children[j]);
			}
		}
	}

	current_level = tmp;
	createCurrentLevelStructure();
	 */
}


void HierarchicalChromosome::createCurrentLevelStructure() {
	chr.clear();
	for (string c: chrs) {
		chr[c].init();
		for (unsigned int i = 0; i < current_level[c].size(); ++i) {
			chr[c].points.push_back(clusters[current_level[c][i]].pos);
			chr[c].genomic_position.push_back(clusters[current_level[c][i]].genomic_pos);
		}
		chr[c].updateSize();
	}
}

void HierarchicalChromosome::smoothSpline(int n) {

	createCurrentLevelStructure();

	for (string c: chrs) {
		std::vector<vector3> v;
		for (unsigned int i = 0; i < chr[c].points.size(); ++i) v.push_back(chr[c].points[i]);

		std::vector<vector3> vi = interpolateSpline(v, n);
		//std::vector<vector3> vi = interpolateSplineCentripetal(v, n);

		chr_smooth[c].init();
		for (unsigned int i = 0; i < vi.size(); ++i) chr_smooth[c].points.push_back(vi[i]);
		chr_smooth[c].updateSize();
	}

	smooth_factor = n;
}

Chromosome HierarchicalChromosome::createEqidistantModel(int resolution_bp, string chr_ind) {
	createCurrentLevelStructure();

	Chromosome ch;

	if (chrs.size() == 0) error("No chromosomes in the model!");
	if (chr_ind == "") chr_ind = chrs[0];

	// safety checks
	if (std::find(chrs.begin(), chrs.end(), chr_ind) == chrs.end()) error(ftext("Chromosome [%s] not found!", chr_ind.c_str()));
	int size = current_level[chr_ind].size();
	if (size == 0) error("Chromosome has size==0");


	if (size == 1) {
		ch.points.push_back(clusters[current_level[chr_ind][0]].pos);
		return ch;
	}

	int pos_st = clusters[current_level[chr_ind][0]].genomic_pos;
	int pos_end = clusters[current_level[chr_ind][size-1]].genomic_pos;

	vector3 pt;
	for (int pos = pos_st; pos<=pos_end; pos+=resolution_bp) {
		pt = find3DSmoothPosition(chr_ind, pos);
		ch.points.push_back(pt);
		ch.genomic_position.push_back(pos);
	}
	ch.updateSize();

	return ch;
}

void HierarchicalChromosome::toFile(string filename) {
	std::stringstream out;
	std::ofstream file;
	file.exceptions(std::ofstream::badbit | std::ofstream::failbit);

	try {
		file.open(filename);
	} catch(const std::ofstream::failure &ex) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	this->toStringStream(out);

	file << out.str();
	file.close();
}

void HierarchicalChromosome::toStringStream(std::stringstream & out) {
	out << (int)clusters.size() << ' ' << (int)arcs.factors.size() << ' ' << -1 << ' ' << chrs.size() << '\n';

	for (unsigned int i = 0; i < arcs.factors.size(); ++i) out << arcs.factors[i] << ' ';
	out << '\n';

	for (unsigned int i = 0; i < clusters.size(); ++i) clusters[i].toStringStream(out);
	
	for (string c: chrs) {
		out << c << ' ' << chr_root[c] << ' ' << arcs.arcs_cnt[c] << '\n';
		for (int i = 0; i < arcs.arcs_cnt[c]; ++i) arcs.arcs[c][i].toStringStream(out);
	}
}

void HierarchicalChromosome::toFilePreviousFormat(string filename) {
	std::stringstream out;
	std::ofstream file;
	file.exceptions(std::ofstream::badbit | std::ofstream::failbit);

	try {
		file.open(filename);
	} catch(const std::ofstream::failure &ex) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	string chr = chrs[0];
	out << clusters.size() << ' ' << arcs.arcs[chr].size() << ' ' << chr_root[chr] << ' ' << arcs.factors.size() << '\n';

	for (unsigned int i = 0; i < arcs.factors.size(); ++i) out << arcs.factors[i] << ' ';
	out << '\n';
	for (unsigned int i = 0; i < clusters.size(); ++i) clusters[i].toStringStreamPreviousFileFormat(out);

	int n = arcs.arcs[chr].size();
	InteractionArc *arc;
	for (int i = 0; i < n; ++i) {
		arc = &arcs.arcs[chr][i];
		out << arc->start << ' ' << arc->end << ' ' << arc->score << ' ';
		out << arc->factor << ' ' << arc->genomic_start << ' ' << arc->genomic_end << " 0 0\n";
	}

	file << out.str();
	file.close();
}

void HierarchicalChromosome::fromFile(string filename) {
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

	// try to decide which format are we dealing with
	// second value now is n_factors (which is small), before it was n_arcs (which is big).
	// UPDATE 161212: the number of arcs can be small for domain models, so this method is not reliable. Now checking both
	// the number of clusters and arcs (if many clusters (>100) but only a few loops (<5) -> previous mode, otherwise the
	// 2nd value is the number of factors)
	int n, n_factors, chrs_cnt, tmp;
	in_stream >> n >> n_factors >> tmp >> chrs_cnt;
	in_stream.seekg(0);

	if (n > 100 && n_factors < 5) this->fromStringStream(in_stream);
	else this->fromStringStreamPreviousFormat(in_stream);
	
	file.close();
}


void HierarchicalChromosome::fromStringStream(std::stringstream & in) {

	printf("read chromosome\n");
	clusters.clear();
	current_level.clear();

	int n, n_factors, chrs_cnt, tmp;
	in >> n >> n_factors >> tmp >> chrs_cnt;

	for (int i = 0; i < n_factors; ++i) {
		std::string str;
		in >> str;
		arcs.factors.push_back(str);
	}

	printf("clusters...\n");
	for (int i=0; i<n; i++) {
		Cluster c;
		c.fromStringStream(in);
		clusters.push_back(c);
	}

	printf("chromosomes...\n");
	int n_arcs, root_ind;
	for (int i = 0; i < chrs_cnt; ++i) {
		std::string str;
		in >> str >> root_ind >> n_arcs;
		
		std::cout << str << "...\n";
		chrs.push_back(str);
		chr_root[str] = root_ind;

		for (int j = 0; j < n_arcs; ++j) {
			InteractionArc arc;
			arc.fromStringStream(in);
			arcs.arcs[str].push_back(arc);

			if (arc.start != -1 && arc.end != -1) {
				vector_insert_unique(clusters[arc.start].arcs, j);
				vector_insert_unique(clusters[arc.end].arcs, j);
			}
		}

		arcs.arcs_cnt[str] = n_arcs;
	}

	printf("done!\n");

	/*
	for (int i = 0; i < clusters.size(); ++i) {
		for (int j = 0; j < clusters[i].children.size(); ++j) {
			clusters[clusters[i].children[j]].parent = i;
		}
	}
	 */
	useTopLevel();
}

void HierarchicalChromosome::fromStringStreamPreviousFormat(std::stringstream & in) {

	printf("old format detected\n");
	clusters.clear();
	current_level.clear();

	int n, n_arcs, n_factors, root;
	in >> n >> n_arcs >> root >> n_factors;

	chrs.push_back("chr");
	chr_root["chr"] = root;

	printf("factors\n");
	for (int i = 0; i < n_factors; ++i) {
		std::string str;
		in >> str;
		arcs.factors.push_back(str);
	}

	printf("clusters\n");
	int st, end, gpos, tmp;
	float x, y, z;
	int children_cnt;

	for (int i=0; i<n; i++) {
		Cluster c;
		//c.fromFile(file);
		in >> gpos >> st >> end >> x >> y >> z >> children_cnt;

		for (int i = 0; i < children_cnt; ++i) {
			in >> tmp;
			c.children.push_back(tmp);
		}

		c.start = st;
		c.end = end;
		c.genomic_pos = gpos;
		c.pos.set(x, y, z);

		clusters.push_back(c);
	}


	//int st, end, score, factor;
	printf("arcs\n");
	int score_raw, factor_raw;

	for (int i = 0; i < n_arcs; ++i) {
		InteractionArc arc;

		in >> arc.start >> arc.end >> arc.score >> arc.factor;
		in >> arc.genomic_start >> arc.genomic_end >> score_raw >> factor_raw;

		arc.eff_score = arc.score;

		arcs.arcs["chr"].push_back(arc);

		if (arc.start != -1 && arc.end != -1) {
			vector_insert_unique(clusters[arc.start].arcs, i);
			vector_insert_unique(clusters[arc.end].arcs, i);
		}
	}
	arcs.arcs_cnt["chr"] = n_arcs;
}

bool HierarchicalChromosome::fromHiCEvo(std::string filename) {
	std::ifstream file;
	file.exceptions(std::ifstream::badbit | std::ifstream::failbit);

	try {
		file.open(filename);
	} catch(const std::ifstream::failure &ex) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return false;
	}

	std::string full_input(std::istreambuf_iterator<char>{file}, {});
	std::stringstream in_stream(full_input);

	clusters.clear();
	current_level.clear();

	arcs.clear();	// clear anchors, arcs, factors, chrs list etc.

	string chr_name = "chr";
	chrs.push_back(chr_name);
	chr_root[chr_name] = 0;

	arcs.arcs_cnt[chr_name] = 0;

	int n, tmp, start, end, par, n_child;
	float x, y, z;

	in_stream >> n;		// # of regions

	for (int i=0; i<n; i++) {
		in_stream >> start >> end >> par >> x >> y >> z >> n_child;

		Cluster c(start, end);
		c.pos.x = x;
		c.pos.y = y;
		c.pos.z = z;
		c.parent = par;

		for (int j = 0; j < n_child; ++j) {
			in_stream >> tmp;
			c.children.push_back(tmp);
		}

		clusters.push_back(c);
	}

	file.close();
	current_level.clear();
	current_level[chr_name].push_back(0);
	return true;
}

bool HierarchicalChromosome::fromTxt(string filename) {
	printf("read chromosome from txt\n");

	std::ifstream file;
	file.exceptions(std::ifstream::badbit | std::ifstream::failbit);

	try {
		file.open(filename);
	} catch(const std::ifstream::failure &ex) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return false;
	}

	std::string full_input(std::istreambuf_iterator<char>{file}, {});
	std::stringstream in_stream(full_input);

	clusters.clear();
	current_level.clear();

	arcs.clear();	// clear anchors, arcs, factors, chrs list etc.

	string chr_name = "chrN";
	chrs.push_back(chr_name);
	chr_root[chr_name] = 0;

	arcs.arcs_cnt[chr_name] = 0;

	int n;
	float x, y, z;

	in_stream >> n;		// # of regions

	Cluster root(1000, 1000*n+100);
	root.children.push_back(1);
	clusters.push_back(root);

	Cluster root_chr(1000, 1000*n+100);
	for (int i=2; i<=n+1; i++) root_chr.children.push_back(i);
	root_chr.parent = 0;
	clusters.push_back(root_chr);

	for (int i=0; i<n; i++) {
		in_stream >> x >> y >> z;
		
		Cluster c(1000*(i+1), 1000*(i+1)+100);
		c.pos.x = x;
		c.pos.y = y;
		c.pos.z = z;
		c.parent = 1;
		clusters.push_back(c);
	}

	file.close();

	current_level.clear();
	current_level[chr_name].push_back(0);
	return true;
}

//vector<int> HierarchicalChromosomeMixed::findFlankingAnchors(int genomic_position) {
//	vector<int> r;
//	for (int i = 0; i < current_level.size(); ++i) {
//		//		if (clusters[current_level[i]].end < start || clusters[current_level[i]].start > end) continue; // outside the range
//		//
//		//		hc.clusters.push_back(clusters[current_level[i]]);
//		//	}
//	}
//	return r;
//}

vector<int> HierarchicalChromosome::findFlankingAnchors(string chr, int genomic_position_start, int genomic_position_end) {
	vector<int> r;
	int p;

	int flank_left = -1, flank_right = -1;
	for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
		//printf("i=%d ", i);
		p = current_level[chr][i];
		//printf("p=%d\n", p);

		if (flank_left == -1) {
			if (clusters[p].end > genomic_position_start) {
				flank_left = p;
				r.push_back(p);
				r.push_back(i);
			}
		}
		else if (flank_right == -1) {
			if (clusters[p].start > genomic_position_end) {
				flank_right = p;
				r.push_back(p);
				r.push_back(i);
				break;
			}
		}
	}
	if (r.size() < 4) r.clear(); // something went wrong
	return r;
}


HierarchicalChromosome HierarchicalChromosome::extractFragment(int start, int end) {

	// TODO: update to genome
	HierarchicalChromosome hc;

	//	useLowestLevel();
	//
	//	for (int i = 0; i < current_level.size(); ++i) {
	//		if (clusters[current_level[i]].end < start || clusters[current_level[i]].start > end) continue; // outside the range
	//
	//		hc.clusters.push_back(clusters[current_level[i]]);
	//	}
	//
	//	Cluster r(hc.clusters[0].start, hc.clusters[hc.clusters.size()-1].end);
	//	for (int i = 0; i < hc.clusters.size(); ++i) {
	//		r.children.push_back(i);
	//		hc.clusters[i].parent = hc.clusters.size();
	//		hc.clusters[i].level = 1;
	//	}
	//
	//
	//	// update arcs
	//	for (int i = 0; i < arcs.raw_arcs.size(); ++i) {
	//
	//		//printf("raw %d %d\n", raw_arcs[i].start, raw_arcs[i].end);
	//		int st = -1, end = -1;
	//		for (int j = 0; j < hc.clusters.size() && (st==-1||end==-1); ++j) {
	//			if (hc.clusters[j].contains(arcs.raw_arcs[i].start)) st = j;
	//			if (hc.clusters[j].contains(arcs.raw_arcs[i].end)) end = j;
	//		}
	//
	//		if (st != -1 && end != -1) {
	//			//printf("new arc: %d %d\n", st, end);
	//			hc.arcs.raw_arcs.push_back(arcs.raw_arcs[i]);
	//
	//			InteractionArc arc(st, end, arcs.raw_arcs[i].score, arcs.raw_arcs[i].factor);
	//			hc.arcs.arcs.push_back(arc);
	//		}
	//	}
	//	hc.arcs.arcs_cnt = hc.arcs.arcs.size();
	//	hc.arcs.factors = arcs.factors;
	//
	//	// add root
	//	hc.root = hc.clusters.size();
	//	hc.clusters.push_back(r);

	return hc;
}

Heatmap HierarchicalChromosome::getDistancesHeatmap() {

	// calc size
	int n = 0;
	for (string chr: chrs) n += current_level[chr].size();

	Heatmap h(n);
	int x = 0, y;
	for (string chr_id: chrs) {
		for (unsigned int i = 0; i < current_level[chr_id].size(); ++i) {
			y = 0;
			for (string chr_id2: chrs) {
				for (unsigned int j = 0; j < current_level[chr_id2].size(); ++j) {
					h.v[x][y] = (clusters[current_level[chr_id][i]].pos - clusters[current_level[chr_id2][j]].pos).length();
					y++;
				}
			}

			x++;
		}
	}
	return h;
}


// calculate some basic stats for chromosomes at different scales, eg. mean & max distance, chromosome diameter etc.
void HierarchicalChromosome::getSpatialDistributionStats() {
	// calculate, at different scales (chromosomal, subchromosomal and cluster):
	// * diameter of chromosomes
	// * distances between chromosomes

	int i,j;
	float d;
	float min, max, avg;
	for (int lvl=0; lvl<=0; lvl++) {
		setLevel(lvl);
		printf("\n* level %d\n", lvl);

		printf("diameters:\n");
		for (string chr_id: chrs) {
			d = chr[chr_id].getDiameter();
			printf("%s: %f\n", chr_id.c_str(), d);
		}

		printf("distance from nucleus center:\n");
		vector3 center = getCenter(true);
		for (string chr_id: chrs) {
			d = 0.0f;
			for (unsigned int i = 0; i < current_level[chr_id].size(); ++i) {
				d += (clusters[current_level[chr_id][i]].pos - center).length();
			}
			d /= current_level[chr_id].size();
			printf("%s: %f\n", chr_id.c_str(), d);
		}

		printf("\ninter chromosomal distances:\n");
		for (string chr_id: chrs) {
			for (string chr_id2: chrs) {
				min = 1e9;
				max = 0.0f;
				avg = 0.0f;

				for (unsigned int i = 0; i < current_level[chr_id].size(); ++i) {
					for (unsigned int j = 0; j < current_level[chr_id2].size(); ++j) {
						d = (clusters[current_level[chr_id][i]].pos - clusters[current_level[chr_id2][j]].pos).length();
						if (d < min) min = d;
						if (d > max) max = d;
						avg += d;
					}
				}

				avg /= (float)(current_level[chr_id].size() * current_level[chr_id2].size());
				printf("%f ", avg);
			}
			printf("\n");
		}
	}
}

// create heatmap with pairwise distances between beads for a given chromosome
// level denotes on which hierarchy level we calculate distances:
//   0 - chromosome (ie. calc distances between chromosomal-beads)
//   1 - segment
//   2 - subanchor
// 'chr' is ignored for level==0, otherwise it denote the chromosome we calculate the distances
// we assume that the structure is at the correct level
Heatmap HierarchicalChromosome::createStructuralHeatmap(string chr, int level) {
	Heatmap h;
	float d = 0.0f;
	if (level == 0) {
		int n = chrs.size();
		h.init(n);
		// i,j - indices of chromosomes
		for (int i = 0; i < n; ++i) {
			for (int j = i+1; j < n; ++j) {
				d = (clusters[current_level[chrs[i]][0]].pos-clusters[current_level[chrs[j]][0]].pos).length();
				h.v[i][j] = h.v[j][i] = d;
			}
		}
	}
	else if (level == 1) {
		int n = current_level[chr].size();
		h.init(n);
		// i,j - indices of segment beads
		for (int i = 0; i < n; ++i) {
			for (int j = i+1; j < n; ++j) {
				d = (clusters[current_level[chr][i]].pos-clusters[current_level[chr][j]].pos).length();
				h.v[i][j] = h.v[j][i] = d;
			}
		}
	}
	else if (level == 2) {
		// only consider subanchor beads within the same segment
		//// first, create a map subanchor_id -> segment_id
		//// then for each pair

		// find the total heatmap size
		int n = 0;
		for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
			n += clusters[current_level[chr][i]].children.size();
		}

		h.init(n);

		// now create the heatmap, segment by segment

		int shift = 0;
		for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
			int m = clusters[current_level[chr][i]].children.size(); // size of the current segment
			// for each segment iterate over pairs of subanchor beads
			for (int j = 0; j < m; ++j) {
				int ind1 = clusters[current_level[chr][i]].children[j];
				for (int k = j+1; k < m; ++k) {
					int ind2 = clusters[current_level[chr][i]].children[k];
					d = (clusters[ind1].pos-clusters[ind2].pos).length();
					h.v[shift+j][shift+k] = h.v[shift+k][shift+j] = d;
				}
			}
			shift += m;
		}
	}
	return h;
}

float HierarchicalChromosome::calcDistance(HierarchicalChromosome hc, int level) {
	float ret = 0.0f;
	if (level == 0) {
		// level==0 means chromosome level.
	}
	else {

		for (string c: chrs) {

		}
	}
	return ret;
}

vector3 HierarchicalChromosome::find3DPosition(string chr_ind, int pos) {
	//printf("find %s %d\n", chr_ind.c_str(), pos);
	if (std::find(chrs.begin(), chrs.end(), chr_ind) == chrs.end()) {
		printf("error: no such chromosome! [%s]\n", chr_ind.c_str());
		return vector3(0.0, 0.0, 0.0);
	}
	int size = current_level[chr_ind].size();
	if (size == 0) {
		printf("error: current level has size 0!\n");
		return vector3(0.0f, 0.0f, 0.0f);
	}

	int p;

	// border cases ('pos' is before the first or after the last bead)
	// we can extend it to check if 'pos' is before the *end* of the first, or after the *start* of the end
	p = current_level[chr_ind][0];
	if (pos < clusters[p].end) return clusters[p].pos;
	p = current_level[chr_ind][size-1];
	if (pos >= clusters[p].start) return clusters[p].pos;

	int prev_p = -1;
	for (int i = 0; i < size; ++i) {
		p = current_level[chr_ind][i]; // index of i-th cluster
		// check if 'pos' is contained within
		if (pos >= clusters[p].start && pos <= clusters[p].end) return clusters[p].pos;

		// if 'pos' is within 2 beads interpolate between them
		if (pos < clusters[p].start) {
			float st = (pos - clusters[prev_p].end) / (float)(clusters[p].start - clusters[prev_p].end);
			return interpolate(clusters[prev_p].pos, clusters[p].pos, st);
		}

		prev_p = p;
	}
	return vector3(0.0f, 0.0f, 0.0f); // It will never happen - just a compilation warninng suppression
}

vector3 HierarchicalChromosome::find3DSmoothPosition(string chr_ind, int pos) {

	//safety checks
	if (std::find(chrs.begin(), chrs.end(), chr_ind) == chrs.end()) return vector3();

	int size = current_level[chr_ind].size();
	if (size == 0) return vector3();

	int p;

	// border cases ('pos' is before the first or after the last bead)
	// we can extend it to check if 'pos' is before the *end* of the first, or after the *start* of the end
	p = current_level[chr_ind][0];
	if (pos < clusters[p].end) return clusters[p].pos;
	p = current_level[chr_ind][size-1];
	if (pos >= clusters[p].start) return clusters[p].pos;

	int prev_p = -1;
	for (int i = 0; i < size; ++i) {
		p = current_level[chr_ind][i]; // index of i-th cluster

		// check if 'pos' is contained within
		//if (pos >= clusters[p].start && pos <= clusters[p].end) return clusters[p].pos;

		// if 'pos' is within 2 beads interpolate between them
		if (pos < clusters[p].genomic_pos) {

			// p1, p2 - points we are currently between which
			// p0, p3 - flanking points
			vector3 p1 = clusters[prev_p].pos;
			vector3 p2 = clusters[p].pos;
			vector3 p0 = i>1 ? clusters[current_level[chr_ind][i-2]].pos : mirrorPoint(p1, p2);
			vector3 p3 = i+1<current_level[chr_ind].size() ? clusters[current_level[chr_ind][i+1]].pos : mirrorPoint(p2, p1);


			float t = (pos - clusters[prev_p].genomic_pos) / (float)(clusters[p].genomic_pos - clusters[prev_p].genomic_pos);

			//printf("%d (%d %d) p=%d t=%f \n", pos, clusters[prev_p].genomic_pos, clusters[p].genomic_pos, p, t);
			vector3 r = interpolateSpline(t, p0, p1, p2, p3);
			return r;

			//float st = (pos - clusters[prev_p].end) / (float)(clusters[p].start - clusters[prev_p].end);
			//return interpolate(clusters[prev_p].pos, clusters[p].pos, st);
		}

		prev_p = p;
	}
	return vector3(); // It will never happen - just a compilation warninng suppression
}
