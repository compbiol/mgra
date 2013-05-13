#include "mpbgraph.h"

std::ofstream outlog("/dev/null");

///////////////////////////////////////////////////////////////
void MBGraph::add_edges(size_t index, const Genome& genome, const std::unordered_set<orf_t>& blocks) {
	std::string first_vertex; //in chromosome
	std::string current_vertex; // is rightmost vertex of the previous block
        std::string prev_chr;

	for(auto iter = genome.cbegin(); iter != genome.cend(); ++iter) {
		if (blocks.find(iter->second) != blocks.end()) { 
			if (iter->first.first == prev_chr) {
				if (genome.get_sign(iter->second) > 0) {
					local_graph[index].insert(current_vertex, iter->second + "t");
				} else {
					local_graph[index].insert(current_vertex, iter->second + "h");
				}
			} else { 
				// new chromosome detected				
				if (!first_vertex.empty() && genome.isCircular(prev_chr)) {
					local_graph[index].insert(first_vertex, current_vertex);
				}

				if (genome.get_sign(iter->second) > 0) { 
					first_vertex = iter->second + "t"; 
				} else { 
					first_vertex = iter->second + "h";
				} 
				prev_chr = iter->first.first;
			} 
		
			if (genome.get_sign(iter->second) > 0) { 
				current_vertex = iter->second + "h"; 
			} else { 
				current_vertex = iter->second + "t";
			} 
		} 
	}

	if (!first_vertex.empty() && genome.isCircular(prev_chr)) {
		local_graph[index].insert(first_vertex, current_vertex);
	}
}

Mcolor MBGraph::add_tree(const std::string& tree, std::vector<std::string>& output) {
	if (tree[0] == '(') {
		//non-trivial tree
		if (tree[tree.size() - 1] != ')') {
			std::cerr << "ERROR: Malformed input (sub)tree 1" << std::endl;
			exit(3);
		}

		int p = 0;
		for(size_t j = 1; j < tree.size() - 1; ++j) {
			if (tree[j] == '(') { 
				++p; 
			} else if (tree[j] == ')') {
				--p;
			} else if (tree[j] == ',') { 
				if (p == 0) { 
					Mcolor Q1 = add_tree(tree.substr(1, j - 1), output);
		    			Mcolor Q2 = add_tree(tree.substr(j + 1, tree.size() - j - 2), output);

					DiColor.insert(Q1);
					DiColor.insert(Q2);

					Mcolor Q(Q1, Q2, Mcolor::Union);
		    
					output.push_back("\t\"" + genome_match::mcolor_to_name(Q) + "\"\t->\t\"" + genome_match::mcolor_to_name(Q1) + "\";");
					output.push_back("\t\"" + genome_match::mcolor_to_name(Q) + "\"\t->\t\"" + genome_match::mcolor_to_name(Q2) + "\";");

					return Q;
				} 
			} 
			if (p < 0) {
				std::cerr << "ERROR: Malformed input (sub)tree 2" << std::endl;
				exit(3);
			}
		}
		if (p != 0) {
			std::cerr << "ERROR: Malformed input (sub)tree 3" << endl;
			exit(3);
		}
	} else {
		//single node
		Mcolor Q;
		for(size_t j = 0; j < tree.size(); ++j) {
			std::string c = tree.substr(j, 1);
			if (!genome_match::member_name(c)) {
				std::cerr << "ERROR: Unknown genome in (sub)tree: " << tree << std::endl;
				exit(3);
			}
			Q.insert(genome_match::get_number(c));
		}
		return Q;
	}
}

MBGraph::MBGraph(const std::vector<Genome>& genomes, const ProblemInstance& cfg) 
:SplitBadColors(false)
{
	build_graph(genomes);
	parsing_tree(genomes, cfg);

	if (!cfg.get_target().empty()) { 
		DiColor.erase(genome_match::name_to_mcolor(cfg.get_target()));
	} 

	//check consistency
	for (auto id = DiColor.cbegin(); id != DiColor.cend(); ++id) {
		for(auto jd = id; jd != DiColor.end(); ++jd) {
			Mcolor C(*id, *jd, Mcolor::Intersection);
			if (!C.empty() && C.size() != id->size() && C.size() != jd->size()) {
				std::clog << "Multicolors " << genome_match::mcolor_to_name(*id) << " " << genome_match::mcolor_to_name(*jd) << " have nontrivial intersection, removing the latter" << std::endl;
				DiColor.erase(jd++);
				--jd;
			}
		}
	}
	    
    
	TColor.resize(DiColor.size());

	std::clog << "vecT-consistent colors: " << DiColor.size() << std::endl;

	size_t col = 1;
	for (auto id = DiColor.begin(); id != DiColor.end(); ++id) {
		std::clog << "\t" << genome_match::mcolor_to_name(*id);
		all_T_color.insert(*id);

		// compute complement to *id
		Mcolor C;
		for(size_t j = 0; j < genomes.size(); ++j) {
			if (!(id->mymember(j))) { //IS HERE MEMBER
				C.insert(j);
			} 
		}
		all_T_color.insert(C);

		TColor[col-1] = *id;

		col++;
	}
	std::clog << std::endl;

	// tell where the root resides
	std::clog << "the root resides in between:";
	auto T = DiColor;
	for (auto it = T.begin(); it != T.end(); ++it) {
		for (auto jt = it; ++jt != T.end(); ) {
			Mcolor C(*it, *jt, Mcolor::Intersection);
			if (C.size() == it->size()) {
				T.erase(it++);
				jt = it;
				continue;
			}
			if (C.size() == jt->size()) {
				T.erase(jt++);
				--jt;
			}
		}
		std::clog << " " << genome_match::mcolor_to_name(*it);
	}
	std::clog << std::endl;
}

/***********************Good**********************************/
/*build Obverse edge, vertex set and add all edges*/
void MBGraph::build_graph(const std::vector<Genome>& genomes) { 
	local_graph.resize(genomes.size());
	std::unordered_set<orf_t> blocks; 

	for(size_t i = 0; i < genomes.size(); ++i) { 
		for(auto it = genomes[i].cbegin(); it != genomes[i].cend(); ++it) {
			if (blocks.count(it->second) == 0) { 
				obverse_edges.insert(it->second + "t", it->second + "h");
				blocks.insert(it->second);
				vertex_set.insert(it->second + "t"); 
				vertex_set.insert(it->second + "h"); 
			} 
		}
	} 
	
#ifdef VERSION2
	std::cerr << "Count unique vertex: " << vertex_set.size() << std::endl;
#endif

	for(size_t i = 0; i < genomes.size(); ++i) {
		add_edges(i, genomes[i], blocks);
	}	
} 

/*
Parse given trees and save legend.dot
*/
void MBGraph::parsing_tree(const std::vector<Genome>& genomes, const ProblemInstance& cfg) { 
	std::vector<std::string> trees = cfg.get_trees();

	// add terminal branches	
	for (size_t j = 0; j < genomes.size(); ++j) {
		DiColor.insert(Mcolor(j));
	}

	std::vector<std::string> output;
	for(auto it = trees.cbegin(); it != trees.cend(); ++it) {
		Mcolor C = add_tree(*it, output);
		if (C.size() < genomes.size()) { 
			DiColor.insert(C); // complete multicolor is excluded
		} 
	}

	std::ofstream flegend("legend.dot");
	flegend << "digraph Legend {" << std::endl;

	flegend << "\tnode [style=filled];" << std::endl;

	for (size_t j = 0; j < genomes.size(); ++j) {
		flegend << "\t\"" << cfg.get_name(j) << "\"\t[fillcolor=" <<  cfg.get_RGBcolor(cfg.get_RGBcoeff() * j)  << "];" << std::endl;
	} 


	for(auto it = output.cbegin(); it != output.cend(); ++it) {
		flegend << *it << std::endl;
	} 
	flegend << "}" << std::endl;
	flegend.close();
} 

/*
Метод возвращает список смежных ребер вершине u, вида (индекс_вершины, цвет ребра) 
*/
mularcs_t MBGraph::get_adjacent_multiedges(const vertex_t& u) const { 
	if (u == Infty) {
		std::cerr << "mularcs ERROR: Infinite input" << std::endl;
		exit(1);
	}

	mularcs_t output;
	for (int i = 0; i < size_graph(); ++i) {
		if (local_graph[i].defined(u)) { 
			std::pair<multi_hashmap::const_iterator, multi_hashmap::const_iterator> iters = local_graph[i].equal_range(u);
			for (auto it = iters.first; it != iters.second; ++it) { 
				if (output.find(it->second) != output.end()) { 
					output.find(it->second)->second.insert(i);
				} else { 
					output.insert(std::make_pair(it->second, Mcolor(i)));	
				} 
			}
		} else { 
			if (output.find(Infty) != output.end()) { 
				output.find(Infty)->second.insert(i);
			} else { 
				output.insert(std::make_pair(Infty, Mcolor(i)));
			} 
		} 
	}
	return output;
} 

/*multimularcs_t MBGraph::get_adjacent_multiedges_v2(const vertex_t& u) const { 
	if (u == Infty) {
		std::cerr << "mularcs ERROR: Infinite input" << std::endl;
		exit(1);
	}

	multimularcs_t output;
	for (int i = 0; i < size_graph(); ++i) {
		if (LG[i].defined(u)) { 
			std::pair<multi_hashmap::const_iterator, multi_hashmap::const_iterator> iters = LG[i].equal_range(u);
			for (auto it = iters.first; it != iters.second; ++it) { 
				if (output.find(it->second) != output.end()) { 
					output.find(it->second)->second.insert(i);
				} else { 
					output.insert(std::make_pair(it->second, Mcolor(i)));	
				} 
			} 
		} else { 
			if (output.find(Infty) != output.end()) { 
				output.find(Infty)->second.insert(i);
			} else { 
				output.insert(std::make_pair(Infty, Mcolor(i)));
			} 
		} 
	}

	return output;
} 
*/
/*
Метод возвращает список смежных ребер вершине u, вида (индекс_вершины, цвет ребра), 
но если splitBadColor = true то еще и режет его 
*/
multimularcs_t MBGraph::get_adjacent_multiedges_with_split(const vertex_t& u) const { 
	mularcs_t edges = get_adjacent_multiedges(u);

	multimularcs_t output; 
	for(auto im = edges.begin(); im != edges.end(); ++im) {
		if (SplitBadColors && !member(DiColor, im->second) && im->second.size() < size_graph()) {
			auto C = split_color(im->second);
			for(auto ic = C.begin(); ic != C.end(); ++ic) {
				output.insert(std::make_pair(im->first, *ic)); 
	    		}
		} else { 
			output.insert(std::make_pair(im->first, im->second)); 
		}
	}

	return output; 	
} 

/*
SplitColor(Q) представляет Q в виде дизъюнктного объединения T-consistent мультицветов, т.е. Q = Q1 U ... U Qm
где каждый Qi является T-consistent и все они попарно не пересекаются. SplitColor(Q) возвращает множество { Q1, Q2, ..., Qm }
(в частности, когда Q является T-consistent, имеем m=1 и Q1=Q).
Теперь, когда SplitBadColors = true, то и ребро (x,y) имеет мультицвет Q, то MBG.mulcols(x) будет содежать вместо (Q,x) пары:
(Q1,y), (Q2,y), ..., (Qm,y)
*/
std::set<Mcolor> MBGraph::split_color(const Mcolor& Q) const {
    std::set<Mcolor> S;

    if (member(DiColor, Q)) {
	S.insert(Q);
	return S;
    }

    if (!SplitBadColors) {
        return S;
    }

    equivalence <size_t> EQ;
    for(auto iq = Q.cbegin(); iq!= Q.cend(); ++iq) { 
	EQ.addrel(iq->first, iq->first);
    } 

    for(auto ic = DiColor.begin(); ic != DiColor.end(); ++ic) {
	Mcolor C(*ic, Q, Mcolor::Intersection);
	if (C.size() >= 2 && C.size() == ic->size() ) {
	    for (auto iq = C.begin(); iq != C.end(); ++iq) {
		EQ.addrel(iq->first, C.begin()->first);
	    }
	}
    }

    EQ.update();
    std::map <size_t, Mcolor > cls; //FIXME because cls - is Mcolor
    EQ.get_eclasses(cls);
    for(auto ic = cls.cbegin(); ic != cls.cend(); ++ic) {
	S.insert(ic->second);
    }
    return S;
}

/*
*/
void MBGraph::update_complement_color(const std::vector<Mcolor>& colors) {
	for(auto it = colors.begin(); it != colors.end(); ++it) { 
		if (!CColorM.defined(*it)) { 
			Mcolor temp; 
			for(size_t j = 0; j < size_graph(); ++j) { 
				if ((!it->mymember(j))) { 
					temp.insert(j);
				} 
			} 
			CColorM.insert(*it, temp);
		} 
	} 
} 

/*
*/
bool MBGraph::AreAdjacentBranches(const Mcolor& A, const Mcolor & B) const {
	if (all_T_color.find(A) == all_T_color.end() || all_T_color.find(B) == all_T_color.end()) { 
		return false;
	} 

	Mcolor Q1; 
	Mcolor Q2;

	if (A.size() >= B.size()) {
		Q1 = A;
		Q2 = B;
	} else {
		Q1 = B;
		Q2 = A;
	}

	Mcolor C(Q1, Q2, Mcolor::Difference);
	
	if (C.size() == Q1.size() - Q2.size() && member(all_T_color, C)) { 		
		return true;
	} 

	Mcolor M(Q1, Q2, Mcolor::Union); 
	if (M.size() == Q1.size() + Q2.size() && member(all_T_color, M)) { 	
		return true;
	} 

	return false;
}
///////////////////////////////////////////////////////////////
const Mcolor& MBGraph::CColor(const Mcolor& S) {
	if (!CColorM.defined(S)) {
		Mcolor T;
		for (size_t j = 0; j < size_graph(); ++j) {
			if (!S.mymember(j)) {  
				T.insert(j);
			} 
		}
		CColorM.insert(S, T);
	}
	return CColorM[S];
}

Mcolor MBGraph::CColorRep(const Mcolor& c) {
	Mcolor Q = CColor(c);
	if (Q.size() > c.size() || (Q.size() == c.size() && Q > c)) 	
		return c;
	return Q;
}


