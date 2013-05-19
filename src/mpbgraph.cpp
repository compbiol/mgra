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

MBGraph::MBGraph(const std::vector<Genome>& genomes, const ProblemInstance& cfg) 
:SplitBadColors(false)
,colors(genomes.size(), cfg)
{
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
Метод возвращает список смежных ребер вершине u, вида (индекс_вершины, цвет ребра) 
*/
Mularcs MBGraph::get_adjacent_multiedges(const vertex_t& u) const { 
	if (u == Infty) {
		std::cerr << "mularcs ERROR: Infinite input" << std::endl;
		exit(1);
	}

	Mularcs output;
	for (int i = 0; i < size_graph(); ++i) {
		if (local_graph[i].defined(u)) { 
			std::pair<multi_hashmap::const_iterator, multi_hashmap::const_iterator> iters = local_graph[i].equal_range(u);
			for (auto it = iters.first; it != iters.second; ++it) { 
				if (output.find(it->second) != output.cend()) { 
					output.find(it->second)->second.insert(i);
				} else { 
					output.insert(it->second, Mcolor(i));	
				} 
			}
		} else { 
			if (output.find(Infty) != output.cend()) { 
				output.find(Infty)->second.insert(i);
			} else { 
				output.insert(Infty, Mcolor(i));
			} 
		} 
	}
	return output;
} 

/*
Метод возвращает список смежных ребер вершине u, вида (индекс_вершины, цвет ребра), 
но если splitBadColor = true то еще и режет его 
*/
multimularcs_t MBGraph::get_adjacent_multiedges_with_split(const vertex_t& u) const { 
	Mularcs edges = get_adjacent_multiedges(u);

	multimularcs_t output; 
	for(auto im = edges.cbegin(); im != edges.cend(); ++im) {
		if (SplitBadColors && !member(colors.DiColor, im->second) && im->second.size() < size_graph()) {
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

    if (member(colors.DiColor, Q)) {
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

    for(auto ic = colors.DiColor.begin(); ic != colors.DiColor.end(); ++ic) {
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

bool MBGraph::AreAdjacentBranches(const Mcolor& A, const Mcolor & B) const {
	if (!colors.is_T_consistent_color(A) || !colors.is_T_consistent_color(B)) { //FIXME
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
	
	if (C.size() == Q1.size() - Q2.size() && colors.is_T_consistent_color(C)) { 		
		return true;
	} 

	Mcolor M(Q1, Q2, Mcolor::Union); 
	if (M.size() == Q1.size() + Q2.size() && colors.is_T_consistent_color(M)) { 	
		return true;
	} 

	return false;
}

