#ifndef BREAKPOINT_GRAPH_HPP
#define BREAKPOINT_GRAPH_HPP

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
  { 
  } 

  /*
   * We can change graph use event operations only (about events see in event folder). 
   * If record false, your actions didn't save in history, but changed graph. BE CAREFUL! (default : true)
   * Implementation see in "apply_func_impl.hpp" file. 
   */
  void apply(twobreak_t const & twobreak, bool record = true); 
  void apply(clone_t const & clone, bool record = true); 
  void apply(tandem_duplication_t const & tandem_duplication, bool record = true);
  void apply(insdel_t const & insdel, bool record = true);

  /*
   * Methods saying about vertex property. 
   * Simple vertex is vertex with degree 2 and alternating multicolors on two multiedge. 
   * Indel vertex is vertex corresponding insertion or deletion block 
   * Duplication vertex is vertex corresponding duplication block
   * Vertex have self loop if block have reverse tandem duplication
   * Implementation see in "property_func_impl.hpp" file. 
   */
  bool is_simple_vertex(vertex_t const & x) const;
  bool is_indel_vertex(vertex_t const & x) const;  
  bool is_duplication_vertex(vertex_t const & x) const;
  bool is_have_self_loop(vertex_t const & x) const; 

  /*
   * Method calculates a twobreak score for create complete edge on (x, y).
   * Description about twobreak score see in journal paper for MGRA 2.0. 
   * This twobreak score is an upper bound on really score. 
   * Implementation see in "property_func_impl.hpp" file. 
   */
  int calculate_cost(vertex_t const & x, vertex_t const & y) const;

  /*
   * Test edge on mobility property for edge (x, y). 
   * Edge is non-mobile if corresponding multicolor is T-consistent color (not vec{T}-consistent color) or 
   * each vec-T-consistent color can not form in neighborhood verteces. 
   * First method tested edge (x, y) with all multicolor. 
   * Second method tested edge (x, y) with specified multicolor. 
   * Return false if edge is non-mobile.
   * Implementation see in "property_func_impl.hpp" file. 
   */
  bool is_mobility_edge(vertex_t const & x, vertex_t const & y) const;
  bool is_mobility_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const;

  /* 
   * Method calculates a mobilty score for VIEWED edge with multicolor COLOR according REMOVED edge. 
   * Mobility score is number of verteces where algorithm can form multicolor COLOR. 
   * Implementation see in "property_func_impl.hpp" file.   
   */
  size_t mobility_score(edge_t const & viewed, mcolor_t const & color, edge_t const & removed) const;

  /* 
   * Method check that is edge contain the T-consistent (not vec{T}-consistent) multicolor.
   * Implementation see in "property_func_impl.hpp" file.   
   */
  bool is_contain_T_consistent_color(vertex_t const & x, vertex_t const & y) const; 

  /* 
   * Method checks that input twobreak decrease vertex score. 
   * Description about vertex score see in journal paper for MGRA 2.0.
   * Implementation see in "property_func_impl.hpp" file.   
   */
  std::pair<size_t, size_t> is_decrease_verteces_score(twobreak_t const & twobreak) const; 

  /*
   * get_all_multicolor_edge_[with_info]
   * Method return multicolor[s] between two verteces. 
   * [with_info] - return set of T-consistent multicolors corresponding round.   
   * Implementation see in getter_func_impl.hpp file
   */
  mcolor_t get_all_multicolor_edge(vertex_t const & x, vertex_t const & y) const;
  std::set<mcolor_t> get_all_multicolor_edge_with_info(vertex_t const & u, vertex_t const & v, bool with_bad_edge = true) const;

  /* 
   * get_all_adjacent_multiedges_[with_info]
   * Methods return all adjacent verteces with multicolor. 
   * [with_info] - return all adjacent verteces with T-consistent multicolors corresponding round.   
   * About returned structure you can see on structures/mularcs.hpp file. 
   * Implementation see in getter_func_impl.hpp file.
   */
  mularcs_t get_all_adjacent_multiedges(vertex_t const & u) const; 
  mularcs_t get_all_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge = true) const;

  /*
   * Methods split our multi breakpoint graph on components [with target COLOR]. 
   * Methods returned components such as class equivalence (Equivalence class 
   * implementation see in "utility/equivalence.hpp" file).
   * If not_drop_complete_edge = true method drop edge with complete multicolor (default true). 
   * Implementation see in getter_func_impl.hpp file.
   */
  utility::equivalence<vertex_t> split_on_components(bool not_drop_complete_edge = true) const;
  utility::equivalence<vertex_t> split_on_components_with_color(mcolor_t const & color) const;

  /*
   * Method return all edges corresponding insertion and deletion events in graph
   */
  partgraph_t get_bad_edges() const;

  /*
   * If you get identity multi breakpoint graph after your algorithm you can check some properties on graph. 
   * Remove all verteces corresponding cloning 
   * and all edges corresponding postponed deletion have complete color 
   */
  bool check_consistency_graph() const; 

  /* 
   *
   *
   */  
  inline void update_number_of_splits(size_t round) {
    assert(0 < round && round < 4);
    if (round == 3) { 
      number_of_splits = graph.count_local_graphs() + 1; 
    } else { 
      number_of_splits = round;  
    } 
  } 
  
  DECLARE_GETTER(bool, canformQoo, canformQoo) 
  DECLARE_SETTER(bool, canformQoo, canformQoo)
  
  // Fixme need to think
  inline bool is_postponed_deletion(vertex_t const & u, vertex_t const & v) const { 
    return this->postponed_deletions.defined(u, v);
  }

  /*
   * Delegators for work with multigraph without color. 
   * Detailed see in file multigraph.hpp
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
   * Delegators for work with multicolors. 
   * Detailed see in graph_colors.hpp file.  
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
  * Delegators for work with action history. 
  * Detailed see in file graph_history.hpp
  */
  typedef typename HistoryGraph<mcolor_t>::twobreak_citer twobreak_citer;
  typedef typename HistoryGraph<mcolor_t>::twobreak_criter twobreak_criter;
  DECLARE_DELEGATE_VOID_METHOD( history, change_history, change_history )
  DECLARE_CONST_ITERATOR( twobreak_citer, history, cbegin_2break_history, cbegin_2break_history ) 
  DECLARE_CONST_ITERATOR( twobreak_citer, history, cend_2break_history, cend_2break_history ) 
  DECLARE_CONST_ITERATOR( twobreak_criter, history, crbegin_2break_history, crbegin_2break_history ) 
  DECLARE_CONST_ITERATOR( twobreak_criter, history, crend_2break_history, crend_2break_history ) 

  bool canformQ(vertex_t const & vertex, mcolor_t const & color) const;

  bool supercanformQ(vertex_t const & vertex, mcolor_t const & color) const;
private:
  mularcs_t get_all_adjacent_tc_multiedges(vertex_t const & u) const; 
  mularcs_t get_all_adjacent_vtc_multiedges(vertex_t const & u) const; 
  
  /*
   * can incident multiedges of x form multicolor COLOR (don't check T-consistent formation)
   * if return false, then COLOR cannot be formed
   * if true - who knows... 
   */
  bool canformQ2(arc_t const & input_edge, mcolor_t const & color) const;

  void check_changed_vertex(vertex_t const & vertex);
  
  size_t mobility_score_relative_vertex(vertex_t const & x, mcolor_t const & color, vertex_t const & y, vertex_t const & about) const;

private:
  MultiGraph graph; 
  ColorsGraph<mcolor_t> multicolors;
  HistoryGraph<mcolor_t> history;

  size_t number_of_splits; 
  bool canformQoo;  // safe choice, at later stages may change to false

  partgraph_t insertions;
  partgraph_t postponed_deletions;   

  std::unordered_map<vertex_t, std::list<size_t> > mother_verteces;
};

#include "graph/apply_func_impl.hpp"
#include "graph/property_func_impl.hpp"
#include "graph/getter_func_impl.hpp"

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::canformQ(vertex_t const & vertex, mcolor_t const & color) const {
  assert(multicolors.is_vec_T_consistent_color(color));
  if (vertex == Infty) {
    return canformQoo;
  }
   
  mularcs_t mularcs = get_all_adjacent_multiedges_with_info(vertex);

  bool canform = true;
  std::set<arc_t> how_many;
  for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
    mcolor_t inter_color(color, arc->second, mcolor_t::Intersection); 
    if (inter_color.size() > 0 && inter_color.size() < arc->second.size()) { 
      canform = false;
    } else if (!inter_color.empty()) { 
      how_many.insert(*arc);
    } 
  }

  if (!canform && how_many.size() == 2) {
    return (canformQ2(*how_many.begin(), how_many.rbegin()->second) || canformQ2(*how_many.rbegin(), how_many.begin()->second));      
  }

  return canform;
}

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::canformQ2(arc_t const & input_edge, mcolor_t const & color) const { 
  assert(multicolors.is_vec_T_consistent_color(color)); 

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
void BreakpointGraph<mcolor_t>::check_changed_vertex(vertex_t const & vertex) {
  std::queue<vertex_t> queue; 
  
  if (this->mother_verteces.count(vertex) != 0) { 
    queue.push(vertex);
  } 

  while (!queue.empty()) { 
    vertex_t v = queue.front(); 
    queue.pop();
    
    if (this->mother_verteces.count(v) != 0) { 
      mularcs_t const & mularcs = this->get_all_adjacent_multiedges(v); 
      std::list<size_t> & clones = mother_verteces.find(v)->second; 
      bool remove_pseudo_infty = false;  
      
      for (auto clone_ind = clones.rbegin(); clone_ind != clones.rend() && !remove_pseudo_infty;) { 
        auto clone = history.get_clone(*clone_ind);
        mcolor_t const & target_color = clone.get_mcolor();
        bool is_not_deleted = true;  

        for (auto arc = mularcs.cbegin(); arc != mularcs.cend() && is_not_deleted; ++arc) {
          mcolor_t inter_color(arc->second, target_color, mcolor_t::Intersection);

          if (inter_color.size() == target_color.size()) {  
            auto const & mother_edge = clone.get_mother_edge();

            if (clone.is_have_pseudo_vertex()) { 
              for (auto const & color : mother_edge.second) {
                graph.erase_edge(color.first, arc->first, mother_edge.first); 
                if (arc->first != Infty) {  
                  graph.add_edge(color.first, arc->first, Infty);  
                }
              } 

              mcolor_t other_color = this->get_all_multicolor_edge(mother_edge.first, Infty); 
              for (auto const & color : other_color) {
                graph.erase_edge(color.first, mother_edge.first, Infty);  
              }
              graph.erase_vertex(mother_edge.first);

              is_not_deleted = false;
              remove_pseudo_infty = true;
              queue.push(arc->first);
            } else {
              is_not_deleted = false;
              clones.erase(std::next(clone_ind).base()); 
            }
          }  
        }

        if (is_not_deleted) { 
          ++clone_ind; 
        }
      }

      if (remove_pseudo_infty || (mother_verteces.count(v) != 0 && mother_verteces[v].empty())) {
        mother_verteces.erase(v); 
      }
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
