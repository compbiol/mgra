#ifndef STAGE_HPP
#define STAGE_HPP

template<class graph_t>
struct Algorithm<graph_t>::Stage {
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::edge_t edge_t;  
  typedef typename graph_t::arc_t arc_t;  

  explicit Stage(std::shared_ptr<graph_t> const & gr) 
  : graph(gr)
  , pseudo_infinity_vertex(0)
  {
  }

  virtual bool do_action() = 0;
  
  virtual std::string get_name() = 0;
  
  virtual ~Stage() 
  {
  }

protected: 
  typedef std::set<arc_t> set_arc_t;

  void split_by_mobile_property(vertex_t const & v, mularcs_t const & mularcs, 
      set_arc_t& mobiles, set_arc_t& non_mobiles) const;

  mcolor_t get_min_addit_color_for_tc(mcolor_t const & color) const;

protected:
  std::shared_ptr<graph_t> graph;
  size_t pseudo_infinity_vertex;  
};

template<class graph_t>
typename graph_t::mcolor_type Algorithm<graph_t>::Stage::get_min_addit_color_for_tc(mcolor_t const & color) const { 
  mcolor_t min_color = this->graph->get_complete_color();
  for (auto col = this->graph->cbegin_T_consistent_color(); col != this->graph->cend_T_consistent_color(); ++col) {
    if (col->includes(color)) {
      mcolor_t diff_color(*col, color, mcolor_t::Difference);
      if (diff_color.size() < min_color.size()) {
        min_color = diff_color;
      } 
    } 
  } 
  return min_color;
}

template<class graph_t>
void Algorithm<graph_t>::Stage::split_by_mobile_property(vertex_t const & v, mularcs_t const & mularcs, 
      set_arc_t& mobiles, set_arc_t& non_mobiles) const {
  for (auto const & arc : mularcs) { 
    if (this->graph->is_mobility_edge(v, arc.second, arc.first)) { 
      mobiles.insert(arc); 
    } else { 
      non_mobiles.insert(arc);
    }
  }   
}


#endif


