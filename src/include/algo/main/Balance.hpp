#ifndef BALANCE_STAGE_HPP
#define BALANCE_STAGE_HPP

namespace algo { 

template<class graph_pack_t>
struct Balance : public algo::AbsStage<graph_pack_t> { 
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using insdel_t = typename graph_pack_t::insdel_t;
  
  explicit Balance(size_t max_round = 1)
  : AbsStage<graph_pack_t>("Balance graph and remove insertions and deletions event.", "balance", max_round)
  {
  }

  bool run(graph_pack_t & graph_pack) override;
    
private:
  DECL_LOGGER("BalanceStage");
}; 

template<class graph_pack_t>
bool Balance<graph_pack_t>::run(graph_pack_t & graph_pack) { 
  size_t number_indel_event = 0; 
    
  std::unordered_set<vertex_t > processed; 
  for (vertex_t const & a1 : graph_pack.graph) {  
    vertex_t const & a2 = graph_pack.graph.get_obverse_vertex(a1);
    mularcs_t const & mularcs = graph_pack.get_all_adjacent_multiedges(a1);

    if (graph_pack.is_indel_vertex(a1) && (processed.count(a1) == 0) && graph_pack.is_indel_vertex(a2) && mularcs.size() != 0)  {
      processed.insert(a1); processed.insert(a2);
      mcolor_t const & indel_color = mularcs.union_multicolors(); 
      mcolor_t const & bar_indel_color = graph_pack.multicolors.get_complement_color(indel_color);
      assert(indel_color == graph_pack.get_all_adjacent_multiedges(a2).union_multicolors());
      graph_pack.apply(insdel_t(a1, a2, bar_indel_color, true));
      ++number_indel_event; 
    } 
  }  
  
  return (number_indel_event != 0); 
}

}

#endif
