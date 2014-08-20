#ifndef BREAKPOINT_GRAPH_HPP
#define BREAKPOINT_GRAPH_HPP

#include "genome_match.h"

#include "graph/multigraph.hpp"
#include "graph/graph_colors.hpp"
#include "graph/graph_history.hpp"

template<class mcolor_t>
struct BreakpointGraph {
  typedef mcolor_t mcolor_type; 
  typedef typename structure::Mularcs<mcolor_t> mularcs_t; 
  
  typedef std::pair<vertex_t, vertex_t> edge_t;
  typedef std::pair<vertex_t, mcolor_t> arc_t; 

  typedef typename MultiGraph::partgraph_t partgraph_t; 
  
  typedef typename HistoryGraph<mcolor_t>::twobreak_t twobreak_t; 
  typedef typename HistoryGraph<mcolor_t>::clone_t clone_t;            
  typedef typename HistoryGraph<mcolor_t>::insdel_t insdel_t;
  typedef typename HistoryGraph<mcolor_t>::tandem_duplication_t tandem_duplication_t;
  typedef typename HistoryGraph<mcolor_t>::transform_t transform_t;
  
  template <class conf_t>
  BreakpointGraph(std::vector<MultiGraph::genome_t> const & genomes, conf_t const & cfg)
  : graph(genomes)
  , multicolors(cfg)
  , number_of_splits(1)
  , canformQoo(true)
  , is_mobile_irregular_edge(true)
  { 
  } 

  /*
   * We can change graph use Event operations only. Folow methods do this. 
   * Implementation see in apply_func_impl.hpp file. 
   */
  void apply(twobreak_t const & twobreak, bool record = true); 
  void apply(clone_t const & clone, bool record = true); 
  void apply(tandem_duplication_t const & tandem_duplication, bool record = true);
  void apply(insdel_t const & insdel, bool record = true);

  bool is_simple_vertex(vertex_t const & x) const;
  bool is_indel_vertex(vertex_t const & x) const;  
  bool is_duplication_vertex(vertex_t const & x) const;
  bool is_have_self_loop(vertex_t const & x) const; 

  /*
   *
   */
  int calculate_cost(vertex_t const & x, vertex_t const & y) const;

  size_t mobility_score(edge_t const & viewed, mcolor_t const & color, edge_t const & removed) const;

  /*
   * 
   */
  bool is_mobility_edge(vertex_t const & x, vertex_t const & y) const;
  bool is_mobility_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const;

  bool is_contain_T_consistent_color(vertex_t const & x, vertex_t const & y) const; 

  std::pair<size_t, size_t> is_decrease_verteces_score(twobreak_t const & twobreak) const; 

  /*
   * get_all_multicolor_edge_[with_info]
   * all - return pseudo and real edge with all multicolor
   * with_info - return split edges corresponding rounds   
   * Implementation see in property_func_impl.hpp
   */
  mcolor_t get_all_multicolor_edge(vertex_t const & x, vertex_t const & y) const;
  std::set<mcolor_t> get_all_multicolor_edge_with_info(vertex_t const & u, vertex_t const & v, bool with_bad_edge = true) const;

  /* 
   *
  */
  mularcs_t get_all_adjacent_multiedges(vertex_t const & u) const; 
  mularcs_t get_all_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge = true) const;

  /*
   *
   */
  utility::equivalence<vertex_t> split_on_components(bool not_drop_complete_edge = true) const;
  utility::equivalence<vertex_t> split_on_components_with_color(mcolor_t const & color) const;

  partgraph_t get_bad_edges() const;
  bool check_consistency_graph() const; 

  inline void update_number_of_splits(size_t ns) {
    if (ns == 3) { 
      number_of_splits = graph.count_local_graphs() + 1; 
    } else { 
      number_of_splits = ns;  
    } 
  } 
  
  DECLARE_GETTER(bool, canformQoo, canformQoo) 
  DECLARE_SETTER(bool, canformQoo, canformQoo)
  
  inline bool is_postponed_deletion(vertex_t const & u, vertex_t const & v) const { 
    return this->postponed_deletions.defined(u, v);
  }

  /*
   * Delegators for multigraph. See in file multigraph.hpp
   */
  DECLARE_DELEGATE_PARAM_CONST_METHOD(size_t, graph, degree_vertex, vertex_t const &, degree_vertex)  
  DECLARE_DELEGATE_PARAM_CONST_METHOD(partgraph_t const &, graph, get_partgraph, size_t, get_partgraph)  
  DECLARE_DELEGATE_PARAM_CONST_METHOD(vertex_t, graph, get_obverse_vertex, vertex_t const &, get_obverse_vertex) 
  DECLARE_DELEGATE_CONST_METHOD( bool, graph, is_identity, is_identity )
  DECLARE_DELEGATE_CONST_METHOD( size_t, graph, size, size )
 
  typedef MultiGraph::citer verteces_citer;
  DECLARE_CONST_ITERATOR( verteces_citer, graph, begin, cbegin )  
  DECLARE_CONST_ITERATOR( verteces_citer, graph, end, cend )
  DECLARE_CONST_ITERATOR( verteces_citer, graph, cbegin, cbegin )  
  DECLARE_CONST_ITERATOR( verteces_citer, graph, cend, cend )

  /*
   * Delegators for multicolors. See in file graph_colors.hpp
   */
  std::set<mcolor_t> split_color(mcolor_t const & color) const { 
    return multicolors.split_color_on_tc_color(color, graph.count_local_graphs() + 1);
  }

  DECLARE_DELEGATE_PARAM_METHOD(mcolor_t const&, multicolors, get_complement_color, mcolor_t const &, get_complement_color)
  DECLARE_DELEGATE_PARAM_CONST_METHOD(bool, multicolors, is_T_consistent_color, mcolor_t const &, is_T_consistent_color)
  DECLARE_DELEGATE_PARAM_CONST_METHOD(bool, multicolors, is_vec_T_consistent_color, mcolor_t const &, is_vec_T_consistent_color)
  DECLARE_DELEGATE_CONST_METHOD(mcolor_t const &, multicolors, get_complete_color, get_complete_color)
  DECLARE_DELEGATE_CONST_METHOD(mcolor_t const &, multicolors, get_root_color, get_root_color)
  DECLARE_DELEGATE_CONST_METHOD( size_t, multicolors, count_vec_T_consitent_color, count_vec_T_consitent_color )

  typedef typename ColorsGraph<mcolor_t>::citer citer_colors; 
  DECLARE_CONST_ITERATOR( citer_colors, multicolors, cbegin_T_consistent_color, cbegin_T_consistent_color )  
  DECLARE_CONST_ITERATOR( citer_colors, multicolors, cend_T_consistent_color, cend_T_consistent_color )
  DECLARE_CONST_ITERATOR( citer_colors, multicolors, cbegin_vec_T_consistent_color, cbegin_vec_T_consistent_color )  
  DECLARE_CONST_ITERATOR( citer_colors, multicolors, cend_vec_T_consistent_color, cend_vec_T_consistent_color )

  /*
  * Delegators for history. See in file graph_history.hpp
  */
  typedef typename HistoryGraph<mcolor_t>::twobreak_citer twobreak_citer;
  typedef typename HistoryGraph<mcolor_t>::twobreak_criter twobreak_criter;
  DECLARE_DELEGATE_VOID_METHOD( history, change_history, change_history )
  DECLARE_CONST_ITERATOR( twobreak_citer, history, cbegin_2break_history, cbegin_2break_history ) 
  DECLARE_CONST_ITERATOR( twobreak_citer, history, cend_2break_history, cend_2break_history ) 
  DECLARE_CONST_ITERATOR( twobreak_criter, history, crbegin_2break_history, crbegin_2break_history ) 
  DECLARE_CONST_ITERATOR( twobreak_criter, history, crend_2break_history, crend_2break_history ) 

private:
  mularcs_t get_all_adjacent_tc_multiedges(vertex_t const & u) const; 
  mularcs_t get_all_adjacent_vtc_multiedges(vertex_t const & u) const; 
  
  //can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
  //if return false, then Q cannot be formed
  //if true - who knows... 
  bool canformQ(arc_t const & input_edge, mcolor_t const & color) const;

  void check_changed_vertex(vertex_t const & vertex);
  
  size_t mobility_score_relative_vertex(vertex_t const & x, mcolor_t const & color, vertex_t const & y, vertex_t const & about) const;

private:
  MultiGraph graph; 
  ColorsGraph<mcolor_t> multicolors;
  HistoryGraph<mcolor_t> history;

  size_t number_of_splits; 
  bool canformQoo;  // safe choice, at later stages may change to false
  bool is_mobile_irregular_edge; 

  partgraph_t insertions;
  partgraph_t postponed_deletions;   

  std::unordered_map<vertex_t, std::list<size_t> > mother_verteces;
};

#include "graph/apply_func_impl.hpp"
#include "graph/property_func_impl.hpp"
#include "graph/getter_func_impl.hpp"

//can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
//if return false, then Q cannot be formed
//if true - who knows... 
template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::canformQ(arc_t const & input_edge, mcolor_t const & color) const { 
  if (input_edge.first == Infty) {
    return canformQoo;
  }
   
  mularcs_t mularcs = get_all_adjacent_multiedges_with_info(input_edge.first);

  bool canform = true;
  for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
    mcolor_t inter_color(color, arc->second, mcolor_t::Intersection); 
    if (inter_color.size() > 0 && inter_color.size() < arc->second.size()) { 
      canform = false;
    }
  }

  return canform;
}

template<class mcolor_t>
void BreakpointGraph<mcolor_t>::check_changed_vertex(vertex_t const & v) {
  if (this->mother_verteces.count(v) != 0) { 
    std::list<size_t> & clones = mother_verteces.find(v)->second; 
    mularcs_t const & mularcs = this->get_all_adjacent_multiedges(v); 

    for (auto clone_ind = clones.rbegin(); clone_ind != clones.rend();) { 
      auto clone = history.get_clone(*clone_ind);
      mcolor_t const & target_color = clone.get_mcolor();

      assert(mularcs.size() != 0);

      bool flag = true; 

      for (auto const & arc : mularcs) {
        mcolor_t inter_color(arc.second, target_color, mcolor_t::Intersection);

        if (inter_color.size() == target_color.size()) {  
          //auto const & central_edge = clone.get_central_arc(); 
          //vertex_t mother = mother_edge.first;
          
          auto const & mother_edge = clone.get_mother_edge();
          vertex_t const & where = arc.first;
          
          if (clone.is_have_pseudo_vertex()) { 
            for (auto const & color : mother_edge.second) {
              graph.erase_edge(color.first, where, mother_edge.first); 
              if (where != Infty) {  
                graph.add_edge(color.first, where, Infty);  
              }
            } 

            //std::cerr << genome_match::mcolor_to_name(other_color) << std::endl;
            mcolor_t other_color = this->get_all_multicolor_edge(mother_edge.first, Infty); 
            for (auto const & color : other_color) {
              graph.erase_edge(color.first, mother_edge.first, Infty);  
            }
            graph.erase_vertex(mother_edge.first);
          }

          //twobreak_t finish_twobreak(central_edge.first, where, central_edge.second, mother, mother_edge.second); 
          //history.save_stop_clone(*clone_ind, finish_twobreak);

          flag = false;
          clones.erase(--(clone_ind.base()));
          break;
        }  
      }

      if (flag) { 
        ++clone_ind; 
      }
    }

    if (mother_verteces[v].empty()) {
      mother_verteces.erase(v); 
    }
  }
}

template<class mcolor_t>
typename MultiGraph::partgraph_t BreakpointGraph<mcolor_t>::get_bad_edges() const { 
  partgraph_t answer; 

  for (auto const & edge: insertions) { 
    answer.insert(edge.first, edge.second);
  } 

  for (auto const & edge: postponed_deletions) {
    answer.insert(edge.first, edge.second);
  }

  return answer;
} 

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::check_consistency_graph() const {
  if (mother_verteces.size() != 0) {
    return false;
  }

  size_t bad_postponed_deletion = 0;
  for (auto const & edge : this->postponed_deletions) {
    vertex_t const & a1 = edge.first; 
    vertex_t const & a2 = edge.second;

    mcolor_t const & color = get_all_multicolor_edge(a1, a2);
    if (color != multicolors.get_complete_color()) {    
      ++bad_postponed_deletion; 
    } 
  }

  return (bad_postponed_deletion == 0);
}

#endif
