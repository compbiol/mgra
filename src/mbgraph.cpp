#include "mbgraph.h"

void MBGraph::add_edges(size_t index, const Genome& genome, const std::unordered_set<orf_t>& blocks) {
	vertex_t first_vertex; //in chromosome
	vertex_t current_vertex; // is rightmost vertex of the previous block
        vertex_t prev_chr;

	for(auto iter = genome.cbegin(); iter != genome.cend(); ++iter) {
		if (blocks.find(iter->second.first) != blocks.end()) { 
			if (iter->first.first == prev_chr) {
				if (iter->second.second > 0) {
					local_graph[index].insert(current_vertex, iter->second.first + "t");
				} else {
					local_graph[index].insert(current_vertex, iter->second.first + "h");
				}
			} else { 
				// new chromosome detected				
				if (!first_vertex.empty() && genome.isCircular(prev_chr)) {
					local_graph[index].insert(first_vertex, current_vertex);
				}

				if (iter->second.second > 0) { 
					first_vertex = iter->second.first + "t"; 
				} else { 
					first_vertex = iter->second.first + "h";
				} 
				prev_chr = iter->first.first;
			} 
		
			if (iter->second.second > 0) { 
				current_vertex = iter->second.first + "h"; 
			} else { 
				current_vertex = iter->second.first + "t";
			} 
		} 
	}

	if (!first_vertex.empty() && genome.isCircular(prev_chr)) {
		local_graph[index].insert(first_vertex, current_vertex);
	}
}
