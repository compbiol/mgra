#ifndef BALANCE_STAGE_HPP
#define BALANCE_STAGE_HPP

template<class graph_t>
struct Algorithm<graph_t>::Balance : public Algorithm<graph_t>::Stage { 
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::insdel_t insdel_t;
  
  explicit Balance(std::shared_ptr<graph_t> const & graph)
  : Stage(graph)
  {
  }

  bool do_action() override;
    
  std::string get_name() override { 
    return "Balance graph and remove insertions and deletions event."; 
  }

private:
  DECL_LOGGER("BalanceStage");
}; 

template<class graph_t>
bool Algorithm<graph_t>::Balance::do_action() { 
  INFO("Start balance stage for remove insertions and deletions event")
  size_t number_indel_event = 0; 
    
  std::unordered_set<vertex_t > processed; 
  for (vertex_t const &a1 : *(this->graph)) {  
    vertex_t const & a2 = this->graph->get_obverse_vertex(a1);
    mularcs_t const & mularcs = this->graph->get_all_adjacent_multiedges(a1);

    if (this->graph->is_indel_vertex(a1) && (processed.count(a1) == 0) && this->graph->is_indel_vertex(a2) && mularcs.size() != 0)  {
      processed.insert({a1, a2}); 
      
      mcolor_t const & indel_color = mularcs.union_multicolors(); 
      mcolor_t const & bar_indel_color = this->graph->get_complement_color(indel_color);
      assert(indel_color == this->graph->get_all_adjacent_multiedges(a2).union_multicolors());
      
      /*std::set<mcolor_t> const & set_split_indel = this->graph->split_color(indel_color); 
      bool i_tc = false;
      size_t c_vtc_indel = 0; 
      for (mcolor_t const & color: set_split_indel) {  
        if (this->graph->is_vec_T_consistent_color(color)) {
          ++c_vtc_indel; 
        } else { 
          i_tc = true; 
        } 
      }  
      
      std::set<mcolor_t> const & set_split_bar_indel = this->graph->split_color(bar_indel_color); 
      bool bi_tc = false;
      size_t c_vtc_bar_indel = 0;        
      for (mcolor_t const & color: set_split_bar_indel) {  
        
        if (this->graph->is_vec_T_consistent_color(color)) {
          ++c_vtc_bar_indel; 
        } else { 
          bi_tc = true; 
        } 
      } 
      

      if (i_tc || (c_vtc_bar_indel == std::min(c_vtc_indel, c_vtc_bar_indel))) { 
        ; //is_insertion = true;
      } else if (bi_tc || (c_vtc_indel == std::min(c_vtc_indel, c_vtc_bar_indel))) { 
        std::cerr << a1 << " " << a2 << std::endl;
      }*/

      bool is_insertion = true;
      /*if (i_tc || (c_vtc_bar_indel == std::min(c_vtc_indel, c_vtc_bar_indel))) { 
        is_insertion = true;
      } else if (bi_tc || (c_vtc_indel == std::min(c_vtc_indel, c_vtc_bar_indel))) { 
        is_insertion = false;
      } else { 
        assert(false);
      } */ 

      this->graph->apply(insdel_t(a1, a2, bar_indel_color, is_insertion));
      ++number_indel_event; 
    } 
  }
  
  return (number_indel_event != 0); 
}

#endif
