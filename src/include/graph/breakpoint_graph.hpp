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
  
  typedef std::pair<vertex_t, vertex_t> graph_edge_t;
  typedef std::pair<vertex_t, mcolor_t> graph_acr_t; 

  typedef typename HistoryGraph<mcolor_t>::twobreak_t twobreak_t; 
  typedef typename HistoryGraph<mcolor_t>::clone_t clone_t;            
  typedef typename HistoryGraph<mcolor_t>::insertion_t insertion_t;
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

  void print() const;
  /*
   * We can change graph use Event operations only. Folow methods do this. 
   * Implementation see in apply_func_impl.hpp file. 
   */
  void apply(twobreak_t const & twobreak, bool record = true); 
  void apply(clone_t const & clone, bool record = true); 
  void apply(tandem_duplication_t const & tandem_duplication, bool record = true);
  void apply(insertion_t const & insdel, bool record = true);

  bool is_simple_vertex(vertex_t const & x) const;
  bool is_indel_vertex(vertex_t const & x) const;  
  bool is_duplication_vertex(vertex_t const & x) const;
  bool is_have_self_loop(vertex_t const & x) const; 

  size_t mobility_score(graph_edge_t const & viewed, mcolor_t const & color, graph_edge_t const & removed) const;

  /* 
   * Very specific function. 
   * --- Parametrs: 
   * Both edges have mobilty score = 1 and this score calculate on mobilty_score function.  
   * Input edges are equal to mobility_score function. (Method manualy tested result).
   * --- Return: if true - then edge is mobile. 
   * Implementation see in property_func_impl.hpp file. 
   */
#ifdef PSEUDO_EDGE
  bool is_result_mobility_edge(graph_edge_t const & viewed, mcolor_t const & new_color, graph_edge_t const & removed) const;
#endif
  
  std::pair<size_t, size_t> is_decrease_verteces_score(twobreak_t const & twobreak) const; 

  bool is_contain_T_consistent_color(vertex_t const & x, vertex_t const & y) const; 
  /*
   *
   */
  int calculate_cost(vertex_t const & x, vertex_t const & y) const;

  bool is_mobility_edge(vertex_t const & x, vertex_t const & y) const;
  bool is_mobility_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const;

#ifdef PSEUDO_EDGE
  bool is_pseudo_edge(vertex_t const & mother, mcolor_t const & color) const; 
  bool is_pseudo_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const;
#endif

  /*
   * get_[all | real]_multicolor_edge_[with_info]
   * all - return pseudo and real edge with all multicolor
   * real - return only real edge with real multicolor
   * with_info - return split edges corresponding rounds   
   * Implementation see in property_func_impl.hpp
   */
  mcolor_t get_all_multicolor_edge(vertex_t const & x, vertex_t const & y) const;

#ifdef PSEUDO_EDGE
  mcolor_t get_real_multicolor_edge(vertex_t const & x, vertex_t const & y) const;  
#endif

  std::set<mcolor_t> get_all_multicolor_edge_with_info(vertex_t const & u, vertex_t const & v, bool with_bad_edge = true) const;

  /* get_[all | real | clones]_adjacent_multiedges_[with_info]
   *
  */
  mularcs_t get_all_adjacent_multiedges(vertex_t const & u) const; 
  mularcs_t get_all_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge = true) const;

#ifdef PSEUDO_EDGE
  mularcs_t get_real_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge = true) const; 
  mularcs_t get_clones_adjacent_multiedges_with_info(vertex_t const & u, mcolor_t const & color, bool with_bad_edge = true) const; 
#endif  

  /*
   * go to private section
   */
  mularcs_t get_all_adjacent_vtc_multiedges(vertex_t const & u) const; 
  
  /*
   *
   */
  utility::equivalence<vertex_t> split_on_components(bool not_drop_complete_edge = true) const;
  utility::equivalence<vertex_t> split_on_components_with_color(mcolor_t const & color) const;

  // This metod not used. Never tested
#ifdef PSEUDO_EDGE  
  std::pair<vertex_t, vertex_t> get_real_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const; 
#endif

  bool check_consistency_graph() const; 
  edges_t get_bad_edges() const;

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
    return multicolors.split_color_on_tc_color(color, number_of_splits);
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

  //can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
  //if return false, then Q cannot be formed
  //if true - who knows... 
  bool canformQ(graph_acr_t const & input_edge, mcolor_t const & color) const;

  void check_changed_vertex(vertex_t const & vertex);
  
  size_t mobility_score_relative_vertex(vertex_t const & x, mcolor_t const & color, vertex_t const & y, vertex_t const & about) const;

#ifdef PSEUDO_EDGE    
  std::set<mcolor_t> split_edge_on_pseudo_edge(vertex_t const & u, mcolor_t const & color, vertex_t const & v) const;
#endif

private:
  MultiGraph graph; 
  ColorsGraph<mcolor_t> multicolors;
  HistoryGraph<mcolor_t> history;

  size_t number_of_splits; 
  bool canformQoo;  // safe choice, at later stages may change to false
  bool is_mobile_irregular_edge; 

  edges_t insertions;
  edges_t postponed_deletions;   

  std::unordered_map<vertex_t, std::list<size_t> > mother_verteces;
};

#include "graph/apply_func_impl.hpp"
#include "graph/property_func_impl.hpp"
#include "graph/getter_func_impl.hpp"

//can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
//if return false, then Q cannot be formed
//if true - who knows... 
template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::canformQ(graph_acr_t const & input_edge, mcolor_t const & color) const { 
  if (input_edge.first == Infty) {
    return canformQoo;
  }
   
  mularcs_t mularcs;
#ifdef PSEUDO_EDGE
  if (is_pseudo_edge(input_edge.first, input_edge.second)) { 
    mularcs = get_clones_adjacent_multiedges_with_info(input_edge.first, input_edge.second); 
  } else { 
    mularcs = get_real_adjacent_multiedges_with_info(input_edge.first); 
  }
#else 
  mularcs = get_all_adjacent_multiedges_with_info(input_edge.first);
#endif

  bool canform = true;
  for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
    mcolor_t inter_color(color, arc->second, mcolor_t::Intersection); 
    if (inter_color.size() > 0 && inter_color.size() < arc->second.size()) { 
      canform = false;
    }
  }

  return canform;
}

#ifdef PSEUDO_EDGE    
template<class mcolor_t>
std::set<mcolor_t> BreakpointGraph<mcolor_t>::split_edge_on_pseudo_edge(vertex_t const & u, mcolor_t const & color, vertex_t const & v) const { 
  std::set<mcolor_t> result({color});

  auto const & split_lambda = [&](vertex_t const & a) -> void { 
    if (this->mother_verteces.count(a) != 0) {

      auto const & clones = this->mother_verteces.find(a)->second;
      for (auto clone_ind = clones.crbegin(); clone_ind != clones.crend(); ++clone_ind) { 
        std::set<mcolor_t> current;
        mcolor_t const & clone_color = history.get_clone(*clone_ind).get_mcolor();
        for (mcolor_t const & edge_color : result) {  
          mcolor_t inter_color(clone_color, edge_color, mcolor_t::Intersection);
          if (inter_color.size() > 0 && inter_color.size() < clone_color.size() && inter_color.size() != edge_color.size()) { 
            current.insert(inter_color);
            current.insert(mcolor_t(edge_color, inter_color, mcolor_t::Difference));
          } else {
            current.insert(edge_color);
          } 
        }
        result = current; 
      } 
    }
  };         

  split_lambda(u);
  split_lambda(v);
  
  return result;
}
#endif

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
edges_t BreakpointGraph<mcolor_t>::get_bad_edges() const { 
  edges_t answer; 

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
      //std::cerr << a1 << " " << a2 << std::endl;
      ++bad_postponed_deletion; 
    } 
  }
  return (bad_postponed_deletion == 0);
}

template<class mcolor_t> 
void BreakpointGraph<mcolor_t>::print() const {  
  std::cerr << "CLONING " << mother_verteces.size() << std::endl;
  for (auto const & mother : this->mother_verteces) {
    std::cerr << mother.first << " ";
    auto const & clones = mother.second; 
    for (auto const & clone : clones) { 
      std::cerr << genome_match::mcolor_to_name(history.get_clone(clone).get_mcolor()) << " ";
    }
    std::cerr << std::endl;
  }
}


#endif
