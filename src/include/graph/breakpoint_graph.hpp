#ifndef BREAKPOINT_GRAPH_HPP
#define BREAKPOINT_GRAPH_HPP

#include "genome_match.h"

#include "graph/multigraph.hpp"
#include "graph/graph_colors.hpp"
#include "graph/graph_history.hpp"

template<class mcolor_t>
struct BreakpointGraph : public MBGraph {
  typedef mcolor_t mcolor_type; 
  typedef typename HistoryGraph<mcolor_t>::twobreak_t twobreak_t; 
  typedef typename HistoryGraph<mcolor_t>::clone_t clone_t;            
  typedef typename HistoryGraph<mcolor_t>::insertion_t insertion_t;
  typedef typename HistoryGraph<mcolor_t>::tandem_duplication_t tandem_duplication_t;
  typedef typename HistoryGraph<mcolor_t>::transform_t transform_t;
  typedef typename structure::Mularcs<mcolor_t> mularcs_t; 
  
  template <class conf_t>
  BreakpointGraph(std::vector<MBGraph::genome_t> const & genomes, conf_t const & cfg)
  : MBGraph(genomes)
  , multicolors(cfg)
  , number_of_splits(1)
  , canformQoo(true)
  , is_mobile_irregular_edge(true)
  { 
  } 

  void print() const; //REMOVE

////

  void apply(twobreak_t const & twobreak, bool record = true); 
  void apply(clone_t const & clone, bool record = true); 
  void apply(tandem_duplication_t const & tandem_duplication, bool record = true);
  void apply(insertion_t const & insdel, bool record = true);

  bool is_simple_vertex(vertex_t const & v) const;
  bool is_indel_vertex(vertex_t const & v) const;  
  bool is_duplication_vertex(vertex_t const & v) const;
  bool is_have_self_loop(vertex_t const & v) const; 

  size_t calculate_cost(vertex_t const & u, vertex_t const & v) const;
  bool is_mobility_edge(vertex_t const & x, vertex_t const & y) const;
  bool is_mobility_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const;

  mcolor_t get_edge_multicolor(vertex_t const & u, vertex_t const & v) const;
  std::set<mcolor_t> get_edge_multicolor_with_info(vertex_t const & u, vertex_t const & v, bool with_bad_edge = true) const;

  mularcs_t get_adjacent_multiedges(vertex_t const & u) const; 
  mularcs_t get_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge = true) const;

  utility::equivalence<vertex_t> split_on_components(bool not_drop_complete_edge = true) const;
  utility::equivalence<vertex_t> split_on_components_with_color(mcolor_t const & color) const;

  bool check_consistency_graph() const; 
  edges_t get_bad_edges() const;

  inline void update_number_of_splits(size_t ns) {
    if (ns == 3) { 
      number_of_splits = this->count_local_graphs() + 1; 
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
   * Delegators for multicolors. See in file graph_colors.hpp
   */
  std::set<mcolor_t> split_color(mcolor_t const & color) const { 
    return multicolors.split_color(color, number_of_splits);
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

  //can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
  //if return false, then Q cannot be formed
  //if true - who knows... 
  bool canformQ(vertex_t const & current, mcolor_t const & color) const;
private: 
  void check_changed_vertex(vertex_t const & vertex);
  void new_change_history(size_t index, vertex_t const & where);

  std::set<mcolor_t> split_edge_on_pseudo_edge(vertex_t const & u, mcolor_t const & color, vertex_t const & v) const;

private:
  ColorsGraph<mcolor_t> multicolors;
  HistoryGraph<mcolor_t> history;

  size_t number_of_splits; 
  bool canformQoo;  // safe choice, at later stages may change to false
  bool is_mobile_irregular_edge; 

  edges_t insertions;
  edges_t postponed_deletions;   

  std::unordered_map<vertex_t, std::set<mcolor_t> > pseudo_edges;
  std::unordered_map<vertex_t, std::set<std::pair<mcolor_t, size_t> > > pseudo_infinity_verteces;
};

#include "graph/apply_func_impl.hpp"
#include "graph/property_func_impl.hpp"
#include "graph/getter_func_impl.hpp"

//can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
//if return false, then Q cannot be formed
//if true - who knows... 
template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::canformQ(vertex_t const & current, mcolor_t const & color) const { 
  if (current == Infty) {
    return canformQoo;
  }
    
  mularcs_t const & mularcs = get_adjacent_multiedges_with_info(current); 

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
std::set<mcolor_t> BreakpointGraph<mcolor_t>::split_edge_on_pseudo_edge(vertex_t const & u, mcolor_t const & color, vertex_t const & v) const { 
  std::set<mcolor_t> result({color});

  auto const & split_lambda = [&](vertex_t const & a) -> void { 
    if (this->pseudo_edges.count(a) != 0) {  
      auto pseudo_colors = this->pseudo_edges.find(a)->second;
      for (auto const & pseudo_color : pseudo_colors) {
        std::set<mcolor_t> current;
        for (mcolor_t const & edge_color : result) {  
          mcolor_t inter_color(pseudo_color, edge_color, mcolor_t::Intersection);
          if (inter_color.size() > 0 && inter_color.size() < pseudo_color.size() && inter_color.size() != edge_color.size()) { 
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

template<class mcolor_t>
void BreakpointGraph<mcolor_t>::check_changed_vertex(vertex_t const & v) {
  if (this->pseudo_infinity_verteces.count(v) != 0) { 
    mularcs_t const & mularcs = this->get_adjacent_multiedges(v); 
    auto target_colors = this->pseudo_infinity_verteces.find(v)->second; 
 
    for (auto arc = mularcs.cbegin(); arc != mularcs.cend(); ++arc) {
      auto const & color = arc->second;

      for (auto const & target : target_colors) { 
        mcolor_t inter_color(color, target.first, mcolor_t::Intersection);
        
        if (inter_color.size() == target.first.size()) {
          this->new_change_history(target.second, arc->first);
          this->pseudo_infinity_verteces[v].erase(target);
          this->pseudo_edges[v].erase(target.first);

          if (this->pseudo_infinity_verteces.find(v)->second.empty()) {
            this->pseudo_infinity_verteces.erase(v);
            this->pseudo_edges.erase(v);
          }
        }
      } 
    } 
  }

  if (this->pseudo_edges.count(v) != 0) { 
    mularcs_t const & mularcs = this->get_adjacent_multiedges(v); 
    auto target_colors = this->pseudo_edges.find(v)->second; 
    for (auto arc = mularcs.cbegin(); arc != mularcs.cend(); ++arc) {
      auto const & color = arc->second;
      for (auto const & target : target_colors) { 
        mcolor_t inter_color(color, target, mcolor_t::Intersection);
        if (inter_color.size() == target.size()) {
          this->pseudo_edges[v].erase(target);
          if (this->pseudo_edges.find(v)->second.empty()) {
            this->pseudo_edges.erase(v);
          }   
        }
      } 
    } 
  }    
}

template<class mcolor_t>
void BreakpointGraph<mcolor_t>::new_change_history(size_t index, vertex_t const & where) {
  clone_t const & clone = history.get_clone(index);
  auto const & mother_edge = clone.get_mother_edge();
  
  assert(clone.is_have_pseudo_vertex());  
  //std::cerr << mother_edge.first << " Where = " << where << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;
  for (auto const & color : mother_edge.second) {
    this->erase_edge(color.first, where, mother_edge.first); 
    if (where != Infty) {  
      this->add_edge(color.first, where, Infty);  
    }
  } 

  mcolor_t other_color = this->get_edge_multicolor(mother_edge.first, Infty); 
  //std::cerr << genome_match::mcolor_to_name(other_color) << std::endl;
  for (auto const & color : other_color) {
    this->erase_edge(color.first, mother_edge.first, Infty);  
  }
  this->erase_vertex(mother_edge.first);
} 

template<class mcolor_t>
edges_t BreakpointGraph<mcolor_t>::get_bad_edges() const { 
  edges_t answer; 

  for (auto const & edge: insertions) { 
    answer.insert(edge.first, edge.second);
  } 

  for (auto const & edge: this->postponed_deletions) {
    answer.insert(edge.first, edge.second);
  }
  return answer;
} 

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::check_consistency_graph() const {
  size_t bad_postponed_deletion = 0;

  for (auto const & edge : this->postponed_deletions) {
    vertex_t const & a1 = edge.first; 
    vertex_t const & a2 = edge.second;

    mcolor_t const & color = get_edge_multicolor(a1, a2);
    if (color != multicolors.get_complete_color()) {    
      //std::cerr << a1 << " " << a2 << std::endl;
      ++bad_postponed_deletion; 
    } 
  }
  return (bad_postponed_deletion == 0);
}


template<class mcolor_t> 
void BreakpointGraph<mcolor_t>::print() const {  
std::cerr << "CLONING " << std::endl;
  for (auto const & clone : this->pseudo_edges) {
    auto const & colors = clone.second;
    std::cerr << clone.first << " ";
    for (auto const & color : colors) { 
      std::cerr << genome_match::mcolor_to_name(color) << " ";
    }
    std::cerr << std::endl;
  }
}


#endif
