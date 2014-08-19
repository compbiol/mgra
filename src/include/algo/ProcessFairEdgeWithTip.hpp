#ifndef PROCESS_FAIR_EDGE_WITH_TIP_HPP
#define PROCESS_FAIR_EDGE_WITH_TIP_HPP

template<class graph_t>
struct Algorithm<graph_t>::ProcessFairEdgeWithTip : public Algorithm<graph_t>::Stage {
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;
  
  explicit ProcessFairEdges(std::shared_ptr<graph_t> const & graph, writer::Wdots<graph_t, ProblemInstance<mcolor_t> > w)
  : Stage(graph) 
  , writer(w)
  {
  }
  
  bool do_action() override;
  
  std::string get_name() override {
    return "Process fair edges.";
  }

private: 
  bool is_good_twobreaks(std::vector<twobreak_t> const & twobreaks) const;

};

template<class graph_t>
bool Algorithm<graph_t>::ProcessFairEdges::do_action() { 
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
  
    for(vertex_t const & x : *this->graph) {  
      mularcs_t const & mularcs = this->graph->get_all_adjacent_multiedges(x);

      if (this->graph->is_duplication_vertex(x) || (mularcs.begin()->second == this->graph->get_complete_color())) {
        continue;
      }  

      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
        vertex_t const & y = im->first; // Q == im->second - color of central edge

        if (y != Infty || this->graph->is_duplication_vertex(y)) {
          continue;
        } 
 
        if (!this->graph->is_mobility_edge(x, y)) { 
          /*Is pretty special four cycle*/
          mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          mularcs_t mularcs_y; 
          mularcs_y.erase(x); 
          
          for (auto const & arc_x : mularcs_x) { 
            if (this->graph->degree_vertex(arc_x.first) == 2) {
              mularcs_t mularcs_p = this->graph->get_all_adjacent_multiedges();  
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