#ifndef INCREASE_NUMBER_COMPONENTS_HPP
#define INCREASE_NUMBER_COMPONENTS_HPP

template<class graph_t>
struct Algorithm<graph_t>::IncreaseNumberComponents : public Algorithm<graph_t>::Stage {
	
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::edge_t edge_t; 

	typedef typename graph_t::twobreak_t twobreak_t;
	typedef typename utility::equivalence<vertex_t> equiv_t;
	typedef typename std::map<vertex_t, std::set<edge_t> > bridges_t;

	explicit IncreaseNumberComponents(std::shared_ptr<graph_t> const & graph)
	: Stage(graph) 
	{
	}

	bool do_action() override;
	
	std::string get_name() override { 
		return "Increase number of components.";
	}

private:
	void get_specific_edges(bridges_t & regular_edges, bridges_t & irregular_edges, 
    equiv_t & connected_components, mcolor_t const & color) const {
    for(vertex_t const &x : *this->graph) {
      mularcs_t const & mularcs = this->graph->get_all_adjacent_multiedges(x); 
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
};

template<class graph_t>
bool Algorithm<graph_t>::IncreaseNumberComponents::do_action() { 
	bool isChanged = false;
	size_t number_rear = 0; // number of rearrangements 

	do {
    number_rear = 0;
    for(auto vtc = this->graph->cbegin_vec_T_consistent_color(); vtc != this->graph->cend_vec_T_consistent_color(); ++vtc) {
      bool repeat = true;
      while(repeat) {
        repeat = false;
        bridges_t regular_edges; // reg. edges between diff. connected components of color Q
        bridges_t irregular_edges; // irreg. edges of color *vtc
        equiv_t connected_components = this->graph->split_on_components_with_color(*vtc);
        get_specific_edges(regular_edges, irregular_edges, connected_components, *vtc); 
	
        std::unordered_set<vertex_t> processed;
        // reg. edges between diff. connected components of color *vtc
        for (auto const & bridges : regular_edges) {
		  		if (processed.count(bridges.first) != 0) { 
		    		continue;
		  		}	

          if (bridges.second.size() == 4 && irregular_edges[bridges.first].size() == 0) {
            edge_t const & p = *(bridges.second.begin());
            edge_t const & q = *(++bridges.second.begin());

            edge_t const & r = *(++(++bridges.second.begin()));
            edge_t const & t = *(++(++(++bridges.second.begin())));

            if (processed.count(connected_components[p.second]) == 0 && processed.count(connected_components[q.second]) == 0
                && processed.count(connected_components[r.second]) == 0 && processed.count(connected_components[t.second]) == 0) { 
              typedef std::pair<twobreak_t, twobreak_t> twobreaks_t;

              twobreaks_t possible_twobreaks1(twobreak_t(p, q, *vtc), twobreak_t(r, t, *vtc));              
              twobreaks_t possible_twobreaks2(twobreak_t(p, r, *vtc), twobreak_t(q, t, *vtc));
              twobreaks_t possible_twobreaks3(twobreak_t(p, t, *vtc), twobreak_t(q, r, *vtc));
              
              auto const & calc_score_lambda = [&] (twobreaks_t const & possible_twobreaks) -> int {
                auto first_scores = this->graph->is_decrease_verteces_score(possible_twobreaks.first);
                auto second_scores = this->graph->is_decrease_verteces_score(possible_twobreaks.second);
                return ((int)(first_scores.first + second_scores.first) - (int)(first_scores.second + second_scores.second));
              };   
               
              std::vector<int> scores({calc_score_lambda(possible_twobreaks1), 
                calc_score_lambda(possible_twobreaks2), calc_score_lambda(possible_twobreaks3)}); 

              bool flag = false; 
              for(int elem : scores) { 
                if (elem != 0) { 
                  flag = true; 
                  break;
                } 
              }

              if (flag) {
                auto const & apply_lambda = [&] (twobreaks_t const & possible_twobreaks) -> void {
                  this->graph->apply(possible_twobreaks.first);
                  this->graph->apply(possible_twobreaks.second);
                  number_rear += 2;
                  repeat = true;      
                };

                processed.insert({connected_components[p.first], connected_components[p.second], connected_components[q.second], connected_components[r.second], connected_components[t.second]});
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
		    		edge_t const & p = *(bridges.second.begin());
		    		edge_t const & q = *(bridges.second.rbegin());

		    		// N.B. we have CC[p.first] == CC[q.first] == ie->first
        		if (processed.count(connected_components[p.second]) == 0 && processed.count(connected_components[q.second]) == 0) { 
        			processed.insert({connected_components[p.first], connected_components[p.second], connected_components[q.second]});
              //std::cerr << p.first << " " << p.second << " " << q.first << " " << q.second << " " << genome_match::mcolor_to_name(*vtc) << std::endl;
		    			this->graph->apply(twobreak_t(p, q, *vtc));
		    			++number_rear;
		    			repeat = true;      
        		} 
				  } 

		  		// connected component with a single external edge and look for irregular edges
		  		if (bridges.second.size() == 1 && !repeat) {
		  			bool found = false;
		    		edge_t const & p = *(bridges.second.begin());
		    		edge_t q;
	    
		    		for(auto ii = irregular_edges[bridges.first].cbegin(); (ii != irregular_edges[bridges.first].cend()) && !found;) {
	      			edge_t const & ireg_edge = *ii; // edge
	      			mcolor_t color(*vtc, this->graph->get_all_multicolor_edge(p.first, ireg_edge.first), mcolor_t::Union);
					
  						if (color.size() > vtc->size() && this->graph->is_T_consistent_color(color)) {
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
		      		 
		      		this->graph->apply(twobreak_t(p, q, *vtc));
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

#endif
