#ifndef INCREASE_NUMBER_COMPONENTS_HPP
#define INCREASE_NUMBER_COMPONENTS_HPP

namespace algo { 

template<class graph_pack_t>
struct IncreaseNumberComponents : public algo::AbsStage<graph_pack_t> {
	
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using edge_t = typename graph_pack_t::edge_t; 

  using twobreak_t = typename graph_pack_t::twobreak_t;
	using equiv_t = typename utility::equivalence<vertex_t>;
	using bridges_t = typename std::map<vertex_t, std::set<edge_t> >;

	explicit IncreaseNumberComponents(size_t max_rounds = 1)
	: AbsStage<graph_pack_t>("Increase number of components", "increase_components", max_rounds) 
	{
	}

	bool run(graph_pack_t & graph_pack) override;
	
private:
	void get_specific_edges(graph_pack_t const & graph_pack, 
    bridges_t & regular_edges, 
    bridges_t & irregular_edges, 
    equiv_t & connected_components, mcolor_t const & color) const {
    for(vertex_t const & x : graph_pack.graph) {
      mularcs_t const & mularcs = graph_pack.get_all_adjacent_multiedges(x); 
      for(auto const & arc : mularcs) {    
        if (arc.second == color) { 
          if (arc.first == Infty) { 
            irregular_edges[connected_components[x]].insert(std::make_pair(x, arc.first)); 
          } else if (!connected_components.isequiv(x, arc.first)) {
            regular_edges[connected_components[x]].insert(std::make_pair(x, arc.first)); 
          } 
        } 
      }
    }
  }

private:
  DECL_LOGGER("IncreaseNumberComponents");
};

template<class graph_pack_t>
bool IncreaseNumberComponents<graph_pack_t>::run(graph_pack_t & graph_pack) {   
	bool isChanged = false;
	size_t number_rear = 0; // number of rearrangements 

	do {
    number_rear = 0;
    for(auto vtc = graph_pack.multicolors.cbegin_vec_T_consistent_color(); vtc != graph_pack.multicolors.cend_vec_T_consistent_color(); ++vtc) {
      bool repeat = true;
      while(repeat) {
        repeat = false;
        bridges_t regular_edges; // reg. edges between diff. connected components of color Q
        bridges_t irregular_edges; // irreg. edges of color *vtc
        equiv_t connected_components = graph_pack.split_on_components_with_color(*vtc);
        get_specific_edges(graph_pack, regular_edges, irregular_edges, connected_components, *vtc); 

        std::unordered_set<vertex_t> processed;
        // reg. edges between diff. connected components of color *vtc
        for (auto const & bridges : regular_edges) {
		  		if (processed.count(bridges.first) != 0) { 
		    		continue;
		  		}	

          if (bridges.second.size() == 4 && irregular_edges[bridges.first].size() == 0) {
            TRACE("Case 1: four bridges.")
            edge_t const & p = *(bridges.second.begin());
            edge_t const & q = *(++bridges.second.begin());

            edge_t const & r = *(++(++bridges.second.begin()));
            edge_t const & t = *(++(++(++bridges.second.begin())));

            if (processed.count(connected_components[p.second]) == 0 && processed.count(connected_components[q.second]) == 0
                && processed.count(connected_components[r.second]) == 0 && processed.count(connected_components[t.second]) == 0) { 
              using twobreaks_t = std::pair<twobreak_t, twobreak_t>;

              twobreaks_t possible_twobreaks1(twobreak_t(p, q, *vtc), twobreak_t(r, t, *vtc));              
              twobreaks_t possible_twobreaks2(twobreak_t(p, r, *vtc), twobreak_t(q, t, *vtc));
              twobreaks_t possible_twobreaks3(twobreak_t(p, t, *vtc), twobreak_t(q, r, *vtc));
              
              auto const & calc_score_lambda = [&] (twobreaks_t const & possible_twobreaks) -> int {
                auto first_scores = graph_pack.is_decrease_verteces_score(possible_twobreaks.first);
                auto second_scores = graph_pack.is_decrease_verteces_score(possible_twobreaks.second);
                return ((int)(first_scores.first + second_scores.first) - (int)(first_scores.second + second_scores.second));
              };   
               
              std::vector<int> scores;
              scores.push_back(calc_score_lambda(possible_twobreaks1));
              scores.push_back(calc_score_lambda(possible_twobreaks2));
              scores.push_back(calc_score_lambda(possible_twobreaks3)); 

              bool flag = false; 
              for(int elem : scores) { 
                if (elem != 0) { 
                  flag = true; 
                  break;
                } 
              }

              if (flag) {
                auto const & apply_lambda = [&] (twobreaks_t const & possible_twobreaks) -> void {
                  graph_pack.apply(possible_twobreaks.first);
                  graph_pack.apply(possible_twobreaks.second);
                  number_rear += 2;
                  repeat = true;      
                };

                processed.insert(connected_components[p.first]); processed.insert(connected_components[p.second]);
                processed.insert(connected_components[q.second]); processed.insert(connected_components[r.second]);
                processed.insert(connected_components[t.second]);
                if (scores.cbegin() == std::max_element(scores.cbegin(), scores.cend())) { 
                  apply_lambda(possible_twobreaks1);
                } else if ((scores.cbegin() + 1) == std::max_element(scores.cbegin(), scores.cend())) { 
                  apply_lambda(possible_twobreaks2);
                } else if ((scores.cbegin() + 2) == std::max_element(scores.cbegin(), scores.cend())) { 
                  apply_lambda(possible_twobreaks3);
                }
              }
            } 
          }

		  		if (bridges.second.size() == 2 && irregular_edges[bridges.first].size() == 0) {
            TRACE("Case 2: two bridges.")
            
		    		edge_t const & p = *(bridges.second.begin());
		    		edge_t const & q = *(bridges.second.rbegin());

		    		// N.B. we have CC[p.first] == CC[q.first] == ie->first
        		if (processed.count(connected_components[p.second]) == 0 && processed.count(connected_components[q.second]) == 0) { 
        			processed.insert(connected_components[p.first]); processed.insert(connected_components[p.second]);
              processed.insert(connected_components[q.second]);
              //std::cerr << p.first << " " << p.second << " " << q.first << " " << q.second << " " << genome_match::mcolor_to_name(*vtc) << std::endl;
		    			graph_pack.apply(twobreak_t(p, q, *vtc));
		    			++number_rear;
		    			repeat = true;      
        		} 
				  } 

		  		// connected component with a single external edge and look for irregular edges
		  		if (bridges.second.size() == 1 && !repeat) {
            TRACE("Case 2: one bridge and start to find irregular edge.")

		  			bool found = false;
		    		edge_t const & p = *(bridges.second.begin());
		    		edge_t q;
	    
		    		for(auto ii = irregular_edges[bridges.first].cbegin(); (ii != irregular_edges[bridges.first].cend()) && !found;) {
	      			edge_t const & ireg_edge = *ii; // edge
	      			mcolor_t color(*vtc, graph_pack.get_all_multicolor_edge(p.first, ireg_edge.first), mcolor_t::Union);
					
  						if (color.size() > vtc->size() && graph_pack.multicolors.is_T_consistent_color(color)) {
  							// let check what would happen with (be-)edge e=(p.first, irreg.first)
  							// if e is enriched to T-consistent color, great!
  							//std::cerr << "perfect edge is found" << std::endl;
                q = ireg_edge;
  							found = true;
  							irregular_edges[bridges.first].erase(ii++);
  			      } else {
  			      	++ii;
  			      }
			    	}

			    	// we did not find good edge, but there are some candidates - take any
			    	if (!found && irregular_edges[bridges.first].size() == 1) {
		      		//std::cerr << "somewhat good but unique edge is found" << std::endl;
		      		q = *(irregular_edges[bridges.first].cbegin());
		      		irregular_edges.erase(bridges.first);
		      		found = true;
			    	}

			    	if (!found && irregular_edges[bridges.first].size() == 0) {
			    		//std::cerr << "no irregular edges, do fission" << std::endl;
			      	q = std::make_pair(Infty, Infty);
			    		found = true;
			    	}
		    
			    	if (found) {	      
		      		processed.insert({connected_components[p.first], connected_components[p.second]});
		      		
		      		if (q.first != Infty) { 
                processed.insert(connected_components[q.first]);
		      		} 
		      		 
		      		graph_pack.apply(twobreak_t(p, q, *vtc));
		      		++number_rear;			   
		      		repeat = true;
		        }
				  } 
        }
      }
    }

    if (number_rear != 0) { 
		  isChanged = true;
    } 
  } while (number_rear > 0); 

  return isChanged; 
} 

} 

#endif
