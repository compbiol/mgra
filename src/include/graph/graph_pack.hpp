#ifndef GRAPH_PACK_HPP
#define GRAPH_PACK_HPP

#include "structures/config_struct.hpp"

#include "graph/multigraph.hpp"
#include "graph/multicolors.hpp"
#include "graph/history.hpp"

template<class mcolor_t>
struct GraphPack {
  typedef mcolor_t mcolor_type; 

  using mularcs_t = typename structure::Mularcs<mcolor_t>; 
  
  using edge_t = std::pair<vertex_t, vertex_t>;
  using arc_t = std::pair<vertex_t, mcolor_t>; 

  using partgraph_t = typename MultiGraph::partgraph_t; 
  using twobreak_t = typename History<mcolor_t>::twobreak_t; 
  using clone_t =  typename History<mcolor_t>::clone_t;            
  using insdel_t = typename History<mcolor_t>::insdel_t;
  using tandem_duplication_t = typename History<mcolor_t>::tandem_duplication_t;
  using transform_t = typename History<mcolor_t>::transform_t;
  
  GraphPack(std::vector<MultiGraph::genome_t> const & genomes)
  : graph(genomes)
  , pseudo_infinity_vertex(0)
  , number_of_splits(1)
  , canformQoo(true)
  { 
  } 

  /*
   * We recomended change graph uses event operations only (about events see in event folder). 
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
   * If you get identity multi breakpoint graph after your algorithm you can check some properties on graph. 
   * Remove all verteces corresponding cloning 
   * and all edges corresponding postponed deletion have complete color 
   * Implementation see in "property_func_impl.hpp" file.   
   */
  bool is_consistency_graph() const; 

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
   * Method is return a series of two-breaks, which needs for create edge (x, y)
   * with multicolor = COLOR 
   * If method cannot create edge, who don't know, what returned
   */
  transform_t take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const;

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
   * Split color on T-consistent colors.
   * Detailed see function split_color_on_tc_color in multicolors.hpp file.  
   */
  std::set<mcolor_t> split_color(mcolor_t const & color) const { 
    return multicolors.split_color_on_tc_color(color, graph.count_local_graphs() + 1);
  }

  /*
   * Recalculate values which contains in statistics structure.
   * Detailed see function calculate in statistics.hpp file.  
   */
  void update_graph_statistics() {
    stats.calculate(*this);
  }

  /* 
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

  inline vertex_t request_pseudo_infinity_vertex() {
    return ("o" + std::to_string(pseudo_infinity_vertex++) + "o");
  }

  bool is_prosthetic_chromosome(vertex_t const & x, vertex_t const & y) const { 
    return prosthetic_chromosomes.defined(x, y);
  }

  DECLARE_GETTER(bool, canformQoo, canformQoo) 
  DECLARE_SETTER(bool, canformQoo, canformQoo)
    
public: 
  MultiGraph graph; 
  Multicolors<mcolor_t> multicolors;
  History<mcolor_t> history;

  struct Statistics;
  Statistics stats;

private:
  mularcs_t get_all_adjacent_tc_multiedges(vertex_t const & u) const; 
  mularcs_t get_all_adjacent_vtc_multiedges(vertex_t const & u) const; 
  
  /*
   * can incident multiedges of x form multicolor COLOR (don't check T-consistent formation)
   * if return false, then COLOR cannot be formed
   * if true - who knows... 
   */
  bool canformQ(vertex_t const & vertex, mcolor_t const & color) const;
  bool supercanformQ(vertex_t const & vertex, mcolor_t const & color) const;
  std::set<std::tuple<vertex_t, vertex_t, mcolor_t> > canformQ2(vertex_t const & vertex, mcolor_t const & required_color) const;
  void check_changed_vertex(vertex_t const & vertex);

  size_t mobility_score_relative_vertex(vertex_t const & x, mcolor_t const & color, vertex_t const & y, vertex_t const & about) const;

private:
  size_t pseudo_infinity_vertex;
  size_t number_of_splits; 
  bool canformQoo;  // safe choice, at later stages may change to false

  partgraph_t prosthetic_chromosomes; 
  partgraph_t postponed_deletions;
  std::unordered_map<vertex_t, std::list<size_t> > mother_verteces;
};

#include "graph/apply_func_impl.hpp"
#include "graph/property_func_impl.hpp"
#include "graph/getter_func_impl.hpp"
#include "graph/statistics.hpp"

template<class mcolor_t>
bool GraphPack<mcolor_t>::supercanformQ(vertex_t const & vertex, mcolor_t const & required_color) const { 
  assert(multicolors.is_vec_T_consistent_color(required_color));
  
  if (vertex == Infty) {
    return canformQoo;
  }

  //std::cerr << "-- Start canformQ function " << vertex << " " << genome_match::mcolor_to_name(required_color) << std::endl;
  mularcs_t mularcs = get_all_adjacent_multiedges_with_info(vertex);
   
  bool canform = true;
  std::set<std::tuple<vertex_t, vertex_t, mcolor_t> > edges;

  for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
    mcolor_t inter_color(required_color, arc->second, mcolor_t::Intersection); 
    if (inter_color.size() > 0 && inter_color.size() < arc->second.size()) { 
      canform = false;
    } else if (!inter_color.empty()) { 
      edges.insert(std::make_tuple(vertex, arc->first, arc->second)); 
    }
  }

  if (!canform) { 
    //std::cerr << "Cannot form. Return false" << std::endl; 
    return canform; 
  }

  //std::cerr << "--After basic see, we can continue " << edges.size() << std::endl;
  canform = false;
  for (auto edge = edges.cbegin(); edge != edges.cend() && !canform; ++edge) { 
    vertex_t x; 
    vertex_t y; 
    mcolor_t color_edge; 
    std::tie(x, y, color_edge) = *edge;
    mcolor_t bar_color_edge = mcolor_t(required_color, color_edge, mcolor_t::Difference);
    //std::cerr << " Start see on edge " << x << " " << y << " " << genome_match::mcolor_to_name(color_edge) << std::endl;
    //std::cerr << " COLOR WE NEED " << genome_match::mcolor_to_name(bar_color_edge) << std::endl;

    std::set<mcolor_t> bar_colors = multicolors.split_color_on_vtc_color(bar_color_edge);
    std::set<mcolor_t> saved_colors; 

    while (!bar_colors.empty()) { 
      std::set<mcolor_t> next_colors; 

      /*std::cerr << " Start iteration with:  "; 
      for (auto q = bar_colors.cbegin(); q != bar_colors.cend(); ++q) { 
        std::cerr << genome_match::mcolor_to_name(*q) << " ";
      } 
      std::cerr << std::endl;*/

      //std::cerr << "Start to see on bar colors " << std::endl;
      for (auto q = bar_colors.cbegin(); q != bar_colors.cend(); ++q) { 
        bool flag = canformQ(x, *q) && canformQ(y, *q); 

        if (flag) {
          //std::cerr << "Save " << genome_match::mcolor_to_name(*q) << std::endl; 
          saved_colors.insert(*q);
        } else { 
          //std::cerr << "Split on next " << genome_match::mcolor_to_name(*q) << std::endl;
          std::set<mcolor_t> another = multicolors.split_color_on_next_vtc_color(*q);

          if (!another.empty()) { 
            //std::cerr << "Yehh splited " << std::endl;
            next_colors.insert(another.begin(), another.end());
          } 
        }
      }
      //std::cerr << "End to see on bar colors " << std::endl;

      bar_colors = next_colors; 
    }

    mcolor_t union_color; 
    for (auto const & color: saved_colors) { 
      union_color = mcolor_t(union_color, color, mcolor_t::Union);
    }
    if (union_color == bar_color_edge) { 
      canform = true; 
    } else { 
      canform = false; 
    }
  } 

  //std::cerr << "Final RETURN " << canform << std::endl;
  return canform;
}

template<class mcolor_t>
bool GraphPack<mcolor_t>::canformQ(vertex_t const & vertex, mcolor_t const & color) const {
  assert(multicolors.is_vec_T_consistent_color(color));
  
  if (vertex == Infty) {
    return canformQoo;
  }
   
  mularcs_t mularcs = get_all_adjacent_multiedges_with_info(vertex);

  bool canform = true;
  for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
    mcolor_t inter_color(color, arc->second, mcolor_t::Intersection); 
    if (inter_color.size() > 0 && inter_color.size() < arc->second.size()) { 
      canform = false;
    } 
  }

  return canform;
}

/*template<class mcolor_t>
std::set<std::tuple<vertex_t, vertex_t, mcolor_t> > GraphPack<mcolor_t>::canformQ2(vertex_t const & vertex, mcolor_t const & required_color) const { 
  assert(multicolors.is_vec_T_consistent_color(required_color)); 

  if (vertex == Infty) {
    return canformQoo;
  }
   
  mularcs_t mularcs = get_all_adjacent_multiedges_with_info(vertex);
   
  bool canform = true;
  std::set<std::tuple<vertex_t, vertex_t, mcolor_t> > edges;

  for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
    mcolor_t inter_color(required_color, arc->second, mcolor_t::Intersection); 
    if (inter_color.size() > 0 && inter_color.size() < arc->second.size()) { 
      canform = false;
    } else if (!inter_color.empty()) { 
      edges.insert(std::make_tuple(vertex, arc->first, arc->second)); 
    }
  }

  if (!canform) { 
    edges.clear(); 
  }

  return edges;
}
*/

template<class mcolor_t>
void GraphPack<mcolor_t>::check_changed_vertex(vertex_t const & vertex) {
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

/*

    auto left_edges = canformQ2(x, bar_color_edge); 
    auto right_edges = canformQ2(y, bar_color_edge); 

    if (left_edges.empty() || right_edges.empty()) { 
      std::cerr << "One of intersect in big color" << std::endl;
      canform = false; 
    } else { 
      std::set<mcolor_t> last_set; 
      for (auto left_edge = left_edges.cbegin(); left_edge != left_edges.cend(); ++left_edge) { 
        for (auto right_edge = right_edges.cbegin(); right_edge != right_edges.cend(); ++right_edge) { 
          if (std::get<2>(*left_edge) == std::get<2>(*right_edge)) { 
            last_set.insert(std::get<2>(*left_edge));
          } else if (std::get<2>(*left_edge).includes*(std::get<2>(*right_edge))) { 
            last_set.insert(std::get<2>(*right_edge));
          } else if (std::get<2>(*right_edge).includes*(std::get<2>(*right_edge))) { 
            last_set.insert(std::get<2>(*right_edge));
          }
        } 
      }
    }
  } 



  if (canform) {  
    canform = false; 

    std::cerr << "--Start advanced test " << canform << std::endl;
    auto get_output_edges_lambda = [&] (vertex_t const & v, mcolor_t const & target_color) -> std::set<mcolor_t> { 
      std::set<mcolor_t> result; 
      mularcs_t const & mularcs_v = get_all_adjacent_multiedges_with_info(v);
      
      for (auto const & arc : mularcs_v) {
        if (target_color.includes(arc.second)) { 
          result.insert(arc.second);  
        } 
      }

      return result;       
    };

    for (auto edge = edges.cbegin(); edge != edges.cend() && !canform; ++edge) { 
      vertex_t x; 
      vertex_t y; 
      mcolor_t color_edge; 
      std::tie(x, y, color_edge) = *edge;
      mcolor_t bar_color_edge = mcolor_t(require_color, color_edge, mcolor_t::Difference);

      std::cerr << "WE HAVE " << x << " " << genome_match::mcolor_to_name(color_edge) << " " << y << std::endl;
      std::cerr << "WE NEED " << genome_match::mcolor_to_name(bar_color_edge) << std::endl;

      std::set<mcolor_t> left_x = get_output_edges_lambda(x, bar_color_edge);
      std::set<mcolor_t> right_y = get_output_edges_lambda(y, bar_color_edge);

      std::cerr << "Left have" << std::endl;
      for (auto const & color : left_x) { 
        std::cerr << genome_match::mcolor_to_name(color) << std::endl;
      }

      std::cerr << "Right have" << std::endl;
      for (auto const & color : right_y) { 
        std::cerr << genome_match::mcolor_to_name(color) << std::endl;
      }

      std::queue<mcolor_t> queue; 
      std::set<mcolor_t> bar_colors = multicolors.split_color_on_vtc_color(bar_color_edge);
      for (auto const & color : bar_colors) {
        queue.push(color); 
      }

      canform = true;
      std::set<mcolor_t> saved_colors; 
      std::cerr << "Start algorithm" << std::endl;
      while(!queue.empty() && canform) { 
        mcolor_t q = queue.front(); 
        queue.pop();
        std::cerr << "See on " << genome_match::mcolor_to_name(q) << std::endl;

        if (left_x.count(q) != 0 && right_y.count(q) != 0) { 
          std::cerr << "First case " << std::endl;
          saved_colors.insert(q);
        } else if (left_x.count(q) == 0 && right_y.count(q) != 0) { 
          canform = canformQ(x, q);
          std::cerr << "Second case " << canform << std::endl;
          if (canform) { 
            saved_colors.insert(q);
          }
        } else if (left_x.count(q) != 0 && right_y.count(q) == 0) { 
          canform = canformQ(y, q);
          std::cerr << "Third case " << canform << std::endl;
          if (canform) { 
            saved_colors.insert(q);
          } 
        } else { 
          std::cerr << "Fourth case " << std::endl;
          canform = canformQ(x, q) && canformQ(y, q); 

          if (canform) {
            saved_colors.insert(q); 
          } else {
            canform = true; 
            std::set<mcolor_t> another = multicolors.split_color_on_next_vtc_color(q);
            std::cerr << "Split " << genome_match::mcolor_to_name(q) << std::endl;
            for (auto const & color : another) { 
              std::cerr << "we have " << genome_match::mcolor_to_name(color) << std::endl; 
              queue.push(color); 
            } 
            assert(!another.empty());
            assert(*another.begin() != q);

            std::cerr << "Finish" << std::endl;
          }
        }
      }

      mcolor_t union_color; 
      for (auto const & color: saved_colors) { 
        union_color = mcolor_t(union_color, color, mcolor_t::Union);
      }

      if (union_color == bar_color_edge) { 
        canform = true; 
      } else { 
        std::cerr << "Potential good situation" << std::endl;
        assert(canform == false);
      }
    }

    if (!canform) { 
      std::cerr << "YEHH we find" << std::endl;
    }
  } 
*/
#endif
