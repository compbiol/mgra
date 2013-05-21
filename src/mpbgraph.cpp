#include "mpbgraph.h"

std::ofstream outlog("/dev/null");

void MBGraph::add_edges(size_t index, const Genome& genome, const std::unordered_set<orf_t>& blocks) {
	vertex_t first_vertex; //in chromosome
	vertex_t current_vertex; // is rightmost vertex of the previous block
        vertex_t prev_chr;

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

/*
Метод возвращает список смежных ребер вершине u, вида (индекс_вершины, цвет ребра) 
*/
Mularcs<Mcolor> MBGraph::get_adjacent_multiedges(const vertex_t& u, const ColorsGraph<Mcolor>& colors, bool split_bad_colors) const { 
	if (u == Infty) {
		std::cerr << "mularcs ERROR: Infinite input" << std::endl;
		exit(1);
	}

	Mularcs<Mcolor> output;
	for (int i = 0; i < size_graph(); ++i) {
		if (local_graph[i].defined(u)) { 
			std::pair<partgraph_t::const_iterator, partgraph_t::const_iterator> iters = local_graph[i].equal_range(u);
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

	if (split_bad_colors) { 
		Mularcs<Mcolor> split; 
		for(auto im = output.cbegin(); im != output.cend(); ++im) {
			if (!member(colors.DiColor, im->second) && im->second.size() < size_graph()) {
				auto C = im->second.split_color(colors, split_bad_colors);
				for(auto ic = C.begin(); ic != C.end(); ++ic) {
					split.insert(im->first, *ic); 
	    			}
			} else { 
				split.insert(im->first, im->second); 
			}
		}
		return split; 
	}
 
	return output;
} 
