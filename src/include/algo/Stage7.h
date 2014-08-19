#ifndef CLONE_STAGE_HPP
#define CLONE_STAGE_HPP

template<class graph_t>
struct Algorithm<graph_t>::ProcessClone : public Algorithm<graph_t>::Stage {
  typedef Stage base;
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::clone_t clone_t;        
  
  explicit ProcessClone(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  , pseudo_infinity_vertex(0)
  {
  }
  
  bool do_action() override;
  
  std::string get_name() override { 
    return "Process clone situation.";
  }

private: 
  size_t pseudo_infinity_vertex;  
};

template<class graph_t>
bool Algorithm<graph_t>::ProcessClone::do_action() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 

    for (vertex_t const & x: *this->graph) {  
      mularcs_t const & mularcs = this->graph->get_all_adjacent_multiedges(x);
        
      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
        vertex_t const & y = im->first; // Q == im->second - color of central edge

        if (y == Infty) {
          continue;
        }

        if (!this->graph->is_mobility_edge(x, y)) {  
          mularcs_t && mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y);
          mularcs_y.erase(x);
          mularcs_t && mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          bool sligshot = (mularcs_y.size() == 1) && (mularcs_x.size() != 1) && this->graph->is_vec_T_consistent_color(mularcs_y.cbegin()->second);
          for (auto arc = mularcs_x.cbegin(); arc != (mularcs_x.cend()) && sligshot; ++arc) {
            sligshot = this->graph->is_vec_T_consistent_color(arc->second);
          }
          sligshot = sligshot && (mularcs_y.cbegin()->second == mularcs_x.union_multicolors());

          if (sligshot) { 
            vertex_t const & mother = mularcs_y.begin()->first;        

            if (mother == Infty) { 
              std::string pseudo_vertex = "o" + std::to_string(pseudo_infinity_vertex) + "o";
              //std::cerr << "We have fake 2-break " << x << " " << y << " " << Infty << " " << pseudo_vertex << " " << genome_match::mcolor_to_name(mularcs_y.cbegin()->second) << std::endl;
              ++pseudo_infinity_vertex;
              clone_t clone(arc_t(x, y), mularcs_x, edge_t(pseudo_vertex, mularcs_y.cbegin()->second), true);  
              this->graph->apply(clone);
              ++number_rear;
              found = true;
              assert(this->graph->get_all_multicolor_edge(x, y) == this->graph->get_complete_color()); 
            } else {
              //std::cerr << "Create clone " << x << " " << y << " " << mother << " " 
              //  << genome_match::mcolor_to_name(mularcs_y.cbegin()->second) << std::endl;
              clone_t clone(arc_t(x, y), mularcs_x, *(mularcs_y.cbegin()), false);  
              this->graph->apply(clone);
              ++number_rear;
              found = true;
              assert(this->graph->get_all_multicolor_edge(x, y) == this->graph->get_complete_color()); 
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
