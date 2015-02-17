#ifndef FOUR_CYCLES_STAGE_HPP
#define FOUR_CYCLES_STAGE_HPP

namespace algo { 

template<class graph_pack_t>
struct ProcessFourCycles : public algo::AbsStage<graph_pack_t> { 

  using mcolor_t = typename graph_pack_t::mcolor_type;
  using edge_t = typename graph_pack_t::edge_t;  
  using arc_t = typename graph_pack_t::arc_t; 
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using twobreak_t = typename graph_pack_t::twobreak_t;
  

	explicit ProcessFourCycles(size_t max_round = 1)
	: AbsStage<graph_pack_t>("Process four cycles", "adequate_four_cycles", round_stage_t, max_round) 
	{
	}

	bool run(graph_pack_t& graph_pack) override;
  
private: 
	bool process_good_four_cycles(graph_pack_t& graph_pack);
	//bool process_generalized_good_four_cycles(graph_pack_t& graph_pack);
};

/*
template<class graph_pack_t>
bool ProcessFourCycles<graph_pack_t>::process_generalized_good_four_cycles(graph_pack_t& graph_pack) { 
	bool isChanged = false; 
	auto medians = graph_pack.multicolors.get_medians_colors(); 
	std::vector<edge_t> ancestor_adjacency;
	std::unordered_set<vertex_t> in_ancestor;

	for(auto const & median_color : medians) { 
		mcolor_t left, right, parent;
		std::tie(left, right, parent) = median_color;
		std::string sleft = cfg::get().mcolor_to_name(left); 
		std::string sright = cfg::get().mcolor_to_name(right); 
		std::string sparent = cfg::get().mcolor_to_name(parent);
		TRACE("Start to work with median " << sleft << " " << sright << " " << sparent) 

		std::array<mcolor_t, 3> colors; 
		colors[0] = left; colors[1] = right; colors[2] = parent;

		for (size_t ind_color = 0; ind_color < 3; ++ind_color) { 
			mcolor_t c0 = colors[ind_color]; mcolor_t c1 = colors[(ind_color + 1) % 3]; mcolor_t c2 = colors[(ind_color + 2) % 3]; 
			mcolor_t c01(c0, c1, mcolor_t::Union);
			std::unordered_set<vertex_t> marked = in_ancestor;

			for (vertex_t const & i : graph_pack.graph) { 
				if (graph_pack.graph.degree_vertex(i) != 2) continue;

				vertex_t j = graph_pack.get_all_adjacent_multiedges(i).get_vertex(c01);

				if (j.empty() || j == Infty || graph_pack.graph.degree_vertex(j) != 2) continue;
				if (marked.count(i) != 0 || marked.count(j) != 0) continue; 
				
				vertex_t i2 = graph_pack.get_all_adjacent_multiedges(i).get_vertex(c2); 
				vertex_t j2 = graph_pack.get_all_adjacent_multiedges(j).get_vertex(c2); 

				if (i2.empty() || j2.empty() || i2 == Infty || j2 == Infty) continue;
				if (graph_pack.graph.degree_vertex(i2) != 3 || graph_pack.graph.degree_vertex(j2) != 3) continue;
				if (!c01.includes(graph_pack.get_all_multicolor_edge(i2, j2))) continue;

				mularcs_t mularcs_i2 = graph_pack.get_all_adjacent_multiedges(i2);
				mularcs_i2.erase(i); mularcs_i2.erase(j2); 
				mularcs_t mularcs_j2 = graph_pack.get_all_adjacent_multiedges(j2);
				mularcs_j2.erase(j); mularcs_j2.erase(i2);

				if (!graph_pack.multicolors.is_vec_T_consistent_color(mularcs_i2.begin()->second) 
						|| !graph_pack.multicolors.is_vec_T_consistent_color(mularcs_j2.begin()->second)) continue;

				vertex_t i1 = mularcs_i2.begin()->first; 
				vertex_t j1 = mularcs_j2.begin()->first; 

				if (i1 == Infty || j1 == Infty) continue;
				
				if (graph_pack.get_all_multicolor_edge(i1, j1) != c0) { 
					//INFO("FIND POSSIBLE GOOD CYCLES")
					//INFO(i << " " << j << " " << i2 << " " << j2 << " " << i1 << " " << j1)
					ancestor_adjacency.push_back(edge_t(i, j));
					ancestor_adjacency.push_back(edge_t(i2, j2));				
					in_ancestor.insert(i); in_ancestor.insert(j); in_ancestor.insert(i2); in_ancestor.insert(j2);
					marked.insert(i); marked.insert(j); marked.insert(i2); marked.insert(j2);
				}
			} 
		}
	}

	for (edge_t const & edge : ancestor_adjacency) { 
		if (graph_pack.get_all_multicolor_edge(edge.first, edge.second) != graph_pack.multicolors.get_complete_color()) { 
			std::string central_color = cfg::get().mcolor_to_name(graph_pack.get_all_multicolor_edge(edge.first, edge.second)); 
			TRACE("Apply twobreak on " << edge.first << " " << edge.second << " " << central_color)
			
			mularcs_t mularcs_first = graph_pack.get_all_adjacent_multiedges(edge.first);
			mularcs_first.erase(edge.second);
			mularcs_t mularcs_second = graph_pack.get_all_adjacent_multiedges(edge.second);
			mularcs_second.erase(edge.first);

			for (arc_t const & arc : mularcs_first) { 
				if (graph_pack.multicolors.is_vec_T_consistent_color(arc.second)) { 
					isChanged = true;
					vertex_t v = mularcs_second.get_vertex(arc.second);
					assert(!v.empty());
					graph_pack.apply(twobreak_t(edge.first, arc.first, edge.second, v, arc.second));
				}
			}
		}
	}

	return isChanged;
}
*/

template<class graph_pack_t>
bool ProcessFourCycles<graph_pack_t>::process_good_four_cycles(graph_pack_t& graph_pack) { 
	TRACE("Process good four cycles")
	bool isChanged = false; 
	auto medians = graph_pack.multicolors.get_medians_colors(); 
	std::vector<edge_t> ancestor_adjacency;
	std::unordered_set<vertex_t> in_ancestor;

	for(auto const & median_color : medians) { 
		mcolor_t left, right, parent;
		std::tie(left, right, parent) = median_color;
		std::string sleft = cfg::get().mcolor_to_name(left); 
		std::string sright = cfg::get().mcolor_to_name(right); 
		std::string sparent = cfg::get().mcolor_to_name(parent);
		TRACE("Start to work with median " << sleft << " " << sright << " " << sparent) 

		std::array<mcolor_t, 3> colors; 
		colors[0] = left; colors[1] = right; colors[2] = parent;

		for (size_t ind_color = 0; ind_color < 3; ++ind_color) { 
			mcolor_t c0 = colors[ind_color]; mcolor_t c1 = colors[(ind_color + 1) % 3]; mcolor_t c2 = colors[(ind_color + 2) % 3]; 
			mcolor_t c01(c0, c1, mcolor_t::Union);
			std::unordered_set<vertex_t> marked = in_ancestor;

			for (vertex_t const & i : graph_pack.graph) { 
				if (graph_pack.graph.degree_vertex(i) != 3) continue;

				vertex_t j = graph_pack.get_all_adjacent_multiedges(i).get_vertex(c0);

				if (j.empty() || j == Infty || graph_pack.graph.degree_vertex(j) != 3) continue;
				if (marked.count(i) != 0 || marked.count(j) != 0) continue; 
				
				vertex_t i1 = graph_pack.get_all_adjacent_multiedges(i).get_vertex(c1); 
				vertex_t i2 = graph_pack.get_all_adjacent_multiedges(i).get_vertex(c2); 
				vertex_t j1 = graph_pack.get_all_adjacent_multiedges(j).get_vertex(c1); 
				vertex_t j2 = graph_pack.get_all_adjacent_multiedges(j).get_vertex(c2); 

				TRACE("See on " << i << " " << j << " " << i1 << " " << j1 << " " << i2 << " " << j2)
				if (!i1.empty() && !j2.empty() && i1 != Infty
						&& graph_pack.graph.degree_vertex(i1) == 3 && graph_pack.graph.degree_vertex(j2) == 3 
						&& marked.count(i1) == 0 && marked.count(j2) == 0
						&& graph_pack.get_all_adjacent_multiedges(i1).get_vertex(c0) == j2) { 
					TRACE("FIRST TYPE")
					//INFO("FIRST TYPE")
					ancestor_adjacency.push_back(edge_t(i, i1));
					ancestor_adjacency.push_back(edge_t(j, j2));
					in_ancestor.insert(i); in_ancestor.insert(i1); in_ancestor.insert(j); in_ancestor.insert(j2);
					marked.insert(i); marked.insert(i1); marked.insert(j); marked.insert(j2);
				} else if (!i2.empty() && !j1.empty() && i2 != Infty
						&& graph_pack.graph.degree_vertex(i2) == 3 && graph_pack.graph.degree_vertex(j1) == 3 
						&& marked.count(i2) == 0 && marked.count(j1) == 0
						&& graph_pack.get_all_adjacent_multiedges(i2).get_vertex(c0) == j1) { 
					TRACE("SECOND TYPE")
					//INFO("2 TYPE")
					ancestor_adjacency.push_back(edge_t(i, i2));
					ancestor_adjacency.push_back(edge_t(j, j1));
					in_ancestor.insert(i); in_ancestor.insert(i2); in_ancestor.insert(j); in_ancestor.insert(j1);
					marked.insert(i); marked.insert(i2); marked.insert(j); marked.insert(j1);
				} else if (!i1.empty() && !j1.empty() && i1 != Infty
						&& graph_pack.graph.degree_vertex(i1) == 3 && graph_pack.graph.degree_vertex(j1) == 3 
						&& marked.count(i1) == 0 && marked.count(j1) == 0
						&& graph_pack.get_all_adjacent_multiedges(i1).get_vertex(c2) == j1) { 
					TRACE("THIRD TYPE")	    
					//INFO("3 TYPE")
					ancestor_adjacency.push_back(edge_t(i, j));
					ancestor_adjacency.push_back(edge_t(i1, j1));				
					in_ancestor.insert(i); in_ancestor.insert(j); in_ancestor.insert(i1); in_ancestor.insert(j1);
					marked.insert(i); marked.insert(j); marked.insert(i1); marked.insert(j1);
				} else if (!i2.empty() && !j2.empty() && i2 != Infty
						&& graph_pack.graph.degree_vertex(i2) == 3 && graph_pack.graph.degree_vertex(j2) == 3 
						&& marked.count(i2) == 0 && marked.count(j2) == 0
						&& graph_pack.get_all_adjacent_multiedges(i2).get_vertex(c1) == j2) { 
					TRACE("FOURTH TYPE")
					//INFO("4 TYPE")
					//INFO("See on " << i << " " << j << " " << i1 << " " << j1 << " " << i2 << " " << j2)
					ancestor_adjacency.push_back(edge_t(i, j));
					ancestor_adjacency.push_back(edge_t(i2, j2));				
					in_ancestor.insert(i); in_ancestor.insert(j); in_ancestor.insert(i2); in_ancestor.insert(j2);
					marked.insert(i); marked.insert(j); marked.insert(i2); marked.insert(j2);
				} else if (!i1.empty() && !j1.empty() && !i2.empty() && !j2.empty() && i1 != Infty
						&& graph_pack.graph.degree_vertex(i1) == 3 && graph_pack.graph.degree_vertex(j1) == 3 
						&& marked.count(i1) == 0 && marked.count(j1) == 0 
						&& graph_pack.get_all_adjacent_multiedges(i1).get_vertex(c0) == j1 && (i2 == j1 || j2 == i1)) { 
					TRACE("FIFTH TYPE")
					//INFO("5 TYPE")
		    	ancestor_adjacency.push_back(edge_t(i, i1));
					ancestor_adjacency.push_back(edge_t(j, j1));				
					in_ancestor.insert(i); in_ancestor.insert(j); in_ancestor.insert(i1); in_ancestor.insert(j1);
					marked.insert(i); marked.insert(j); marked.insert(i1); marked.insert(j1);
				} else if (!i1.empty() && !j1.empty() && !i2.empty() && !j2.empty() && i2 != Infty
						&& graph_pack.graph.degree_vertex(i2) == 3 && graph_pack.graph.degree_vertex(j2) == 3 
						&& marked.count(i2) == 0 && marked.count(j2) == 0
						&& graph_pack.get_all_adjacent_multiedges(i2).get_vertex(c0) == j2 && (i1 == j2 || j1 == i2)) { 
					TRACE("SIXTH TYPE")
				 	ancestor_adjacency.push_back(edge_t(i, i2));
					ancestor_adjacency.push_back(edge_t(j, j2));				
					in_ancestor.insert(i); in_ancestor.insert(i2); in_ancestor.insert(j); in_ancestor.insert(j2);
					marked.insert(i); marked.insert(i2); marked.insert(j); marked.insert(j2);
				} else if (!i1.empty() && !j1.empty() && !i2.empty() && !j2.empty() && i1 != Infty && i2 != Infty
						&& graph_pack.graph.degree_vertex(i1) == 3 && graph_pack.graph.degree_vertex(j1) == 3 
						&& graph_pack.graph.degree_vertex(i2) == 3 && graph_pack.graph.degree_vertex(j2) == 3 
						&& marked.count(i1) == 0 && marked.count(j1) == 0 
						&& marked.count(i2) == 0 && marked.count(j2) == 0
						&& graph_pack.get_all_adjacent_multiedges(i1).get_vertex(c0) == j1 
						&& graph_pack.get_all_adjacent_multiedges(i2).get_vertex(c0) == j2 
						&& i1 != j2 && j1 != i2) { 
					TRACE("SEVENTH TYPE")
					//INFO("SEVENTH TYPE")
					ancestor_adjacency.push_back(edge_t(i, j));
					ancestor_adjacency.push_back(edge_t(i1, j1));
					ancestor_adjacency.push_back(edge_t(i2, j2));				
					in_ancestor.insert(i); in_ancestor.insert(j); in_ancestor.insert(i1); in_ancestor.insert(j1);
					in_ancestor.insert(i2); in_ancestor.insert(j2);
					marked.insert(i); marked.insert(j); marked.insert(i1); marked.insert(j1);
					marked.insert(i2); marked.insert(j2);
				} 	
			}
		}		
	}

	for (edge_t const & edge : ancestor_adjacency) { 
		if (graph_pack.get_all_multicolor_edge(edge.first, edge.second) != graph_pack.multicolors.get_complete_color()) { 
			std::string central_color = cfg::get().mcolor_to_name(graph_pack.get_all_multicolor_edge(edge.first, edge.second)); 
			TRACE("Apply twobreak on " << edge.first << " " << edge.second << " " << central_color)
			
			mularcs_t mularcs_first = graph_pack.get_all_adjacent_multiedges(edge.first);
			mularcs_first.erase(edge.second);
			mularcs_t mularcs_second = graph_pack.get_all_adjacent_multiedges(edge.second);
			mularcs_second.erase(edge.first);

			for (arc_t const & arc : mularcs_first) { 
				if (graph_pack.multicolors.is_vec_T_consistent_color(arc.second)) { 
					isChanged = true;
					vertex_t v = mularcs_second.get_vertex(arc.second);
					assert(!v.empty());
					graph_pack.apply(twobreak_t(edge.first, arc.first, edge.second, v, arc.second));
				}
			}
		}
	}

	return isChanged;
}

template<class graph_pack_t>
bool ProcessFourCycles<graph_pack_t>::run(graph_pack_t& graph_pack) {
	bool isChanged = false; 
	bool temp = true;

	while(temp) { 
  	temp = process_good_four_cycles(graph_pack);

  	/*if (!temp) { 
  		temp = process_generalized_good_four_cycles();
  	} */

  	if (temp) { 
  		isChanged = true;
  	} 
	}

	TRACE("Finish process special vec{T}-consistent subgraphs " << isChanged)
  return isChanged;
}

/*
for (vertex_t const & v : *this->graph) { 
			if (graph_pack.is_duplication_vertex(v) || graph_pack.graph.degree_vertex(v) != 3) continue;
		
			std::set<std::vector<vertex_t> > all_possible_four_cycles;
			
			std::function<void(vertex_t const &, vertex_t const &, std::vector<vertex_t> &)> small_dfs = [&] (vertex_t const & x, vertex_t const & prev, std::vector<vertex_t> & cycle) -> void { 
				mularcs_t mularcs_x = graph_pack.get_all_adjacent_multiedges(x); 
				mularcs_x.erase(prev);

				for (arc_t const & arc_x : mularcs_x) { 
					if (graph_pack.is_duplication_vertex(arc_x.first) || graph_pack.graph.degree_vertex(arc_x.first) != 3) continue;
					
					cycle.push_back(arc_x.first);
					if (cycle.size() == 5 && *cycle.begin() == *cycle.rbegin()) { 
						all_possible_four_cycles.insert(cycle);
					} else if (cycle.size() < 5 && arc_x.first != Infty) {  
						small_dfs(arc_x.first, x, cycle);
					} 
					cycle.pop_back();
				}
			};

			std::vector<vertex_t> temp; temp.push_back(v);
			small_dfs(v, "", temp);
			
			if (v == "458h") { 
				std::cerr << all_possible_four_cycles.size() << std::endl;
				for(auto const & cycle : all_possible_four_cycles) { 
					mcolor_t first = graph_pack.get_all_multicolor_edge(cycle[0], cycle[1]);
					mcolor_t second = graph_pack.get_all_multicolor_edge(cycle[1], cycle[3]);
				
					std::stringstream ss; 
					ss << "\nwe have possible good cycle " << cycle[0] << "--" << cycle[1] << "--" << 
						cycle[2] << "--" << cycle[3] << "\n" << cfg::get().mcolor_to_name(first) << " " << cfg::get().mcolor_to_name(second) << "\n";
					std::string temp = ss.str();      					 
					INFO(temp)      					 
				} 
			}

			for(auto const & cycle : all_possible_four_cycles) { 
				mcolor_t first = graph_pack.get_all_multicolor_edge(cycle[0], cycle[1]);
				mcolor_t second = graph_pack.get_all_multicolor_edge(cycle[1], cycle[2]);
				mcolor_t third = graph_pack.get_all_multicolor_edge(cycle[2], cycle[3]);
				mcolor_t four = graph_pack.get_all_multicolor_edge(cycle[3], cycle[4]);

				if (first == third && second == four)  { 
					if ((first == left && second == third) || (first == right && second == left)) {
						mularcs_t mularcs_c0 = graph_pack.get_all_multicolor_edge(cycle[0]);
						mularcs_c0.erase(cycle[3]); mularcs_c0.erase(cycle[1]);
						mularcs_t mularcs_c1 = graph_pack.get_all_multicolor_edge(cycle[1]);
						mularcs_c1.erase(cycle[0]); mularcs_c1.erase(cycle[2]);
						mularcs_t mularcs_c2 = graph_pack.get_all_multicolor_edge(cycle[2]);
						mularcs_c2.erase(cycle[1]); mularcs_c2.erase(cycle[3]);
						mularcs_t mularcs_c3 = graph_pack.get_all_multicolor_edge(cycle[3]);
						mularcs_c3.erase(cycle[2]); mularcs_c3.erase(cycle[0]);

						twobreak_t possible_break11();
						twobreak_t possible_break12();

						twobreak_t possible_break21();
						twobreak_t possible_break22();
						
						/{	 
							std::stringstream ss; 
      				ss << "\nwe have possible good cycle " << cycle[0] << "--" << cycle[1] << "--" << 
      						cycle[2] << "--" << cycle[3] << "\n" << cfg::get().mcolor_to_name(first) << " " << cfg::get().mcolor_to_name(second) << "\n";
							std::string temp = ss.str();      					 
							INFO(temp)      					 
      			}/
					}
				}
			} 
			  
		}
*/

}

#endif 