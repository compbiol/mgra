#ifndef INCREASE_NUMBER_COMPONENTS_HPP
#define INCREASE_NUMBER_COMPONENTS_HPP

template<class graph_t>
struct Algorithm<graph_t>::IncreaseNumberComponents : public Algorithm<graph_t>::Stage {
	typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
	typedef typename graph_t::twobreak_t twobreak_t;
	typedef typename utility::equivalence<vertex_t> equiv_t;
	typedef typename std::map<vertex_t, std::set<arc_t> > bridges_t;

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
      mularcs_t const & mularcs = this->graph->get_adjacent_multiedges(x); 
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

		  		if (bridges.second.size() == 2 && irregular_edges[bridges.first].size() == 0) {
		    		arc_t const & p = *(bridges.second.begin());
		    		arc_t const & q = *(bridges.second.rbegin());

		    		// N.B. we have CC[p.first] == CC[q.first] == ie->first
        		if (processed.count(connected_components[p.second]) == 0 && processed.count(connected_components[q.second]) == 0) { 
        			processed.insert({connected_components[p.first], connected_components[p.second], connected_components[q.second]});
		    			this->graph->apply(twobreak_t(p, q, *vtc));
		    			++number_rear;
		    			repeat = true;      
        		} 
				  } 

		  		// connected component with a single external edge and look for irregular edges
		  		if (bridges.second.size() == 1 && !repeat) {
		  			bool found = false;
		    		arc_t const & p = *(bridges.second.begin());
		    		arc_t q;
	    
		    		for(auto ii = irregular_edges[bridges.first].cbegin(); (ii != irregular_edges[bridges.first].cend()) && !found;) {
	      			arc_t const & ireg_edge = *ii; // edge
	      			mcolor_t color(*vtc, this->graph->get_edge_multicolor(p.first, ireg_edge.first), mcolor_t::Union);
					
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
