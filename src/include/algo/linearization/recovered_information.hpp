#ifndef RECOVERED_INFORMATION_HPP
#define RECOVERED_INFORMATION_HPP

#include "algo/linearization/utils.hpp"
#include "algo/linearization/general_linearize.hpp"

namespace algo { 

template<class graph_pack_t>
struct RecoveredInformation {
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using twobreak_t = typename graph_pack_t::twobreak_t; 
  using transform_t = typename graph_pack_t::transform_t;
  using partgraph_t = typename graph_pack_t::partgraph_t; 

  using genome_t = structure::Genome;
  using chromosome_t = structure::Chromosome; 
  using history_t = std::tuple<transform_t, transform_t, transform_t>;

  struct AncestorInformation { 
    AncestorInformation() = default;

    std::vector<genome_t> genomes;
    std::map<std::pair<mcolor_t, mcolor_t>, transform_t> transformations;
  };

  explicit RecoveredInformation(graph_pack_t const & gp)
  : graph_pack(gp)
  , linearize_algo(gp)
  {
    for (auto const & tree : cfg::get().phylotrees) {
      init_parent_colors(tree.get_root()); 
      if (tree.is_phylogenetic_root()) { 
        mcolor_t const & left_color = tree.get_root()->get_left_child()->get_data();
        mcolor_t const & right_color = tree.get_root()->get_right_child()->get_data();  
        parent_colors.insert(std::make_pair(left_color, right_color));
        parent_colors.insert(std::make_pair(right_color, left_color));
      } 
    }
  }

  /**
   *
   */
  void init_raw_results();

  /**
   *
   */
  void init_linearize_results();

  /**
   *
   */
  void init_target_results();
  
  /**
   *
   */
  AncestorInformation get_results() const { 
    AncestorInformation ret; 

    for (auto const & local_graph : graphs) {
      ret.genomes.push_back(linearize::get_genome(graph_pack, cfg::get().mcolor_to_name(local_graph.first), local_graph.second)); 
    }
    
    for (auto const & elem : transformations) { 
      assert(parent_colors.find(elem.first) != parent_colors.cend());
      ret.transformations.insert(std::make_pair(std::make_pair(elem.first, parent_colors.find(elem.first)->second), elem.second));  
    }
    
    return ret;
  }

private: 
  using tree_t = typename structure::BinaryTree<mcolor_t>; 
  using node_t = typename tree_t::Node; 
  void walk_and_linearize(std::unique_ptr<node_t> const & current);
  void init_parent_colors(std::unique_ptr<node_t> const & current); 

private: 
  graph_pack_t const & graph_pack;
  linearize::GeneralLinearize<graph_pack_t> linearize_algo;

  std::map<mcolor_t, partgraph_t> graphs;
  std::map<mcolor_t, mcolor_t> parent_colors; 
  std::map<mcolor_t, transform_t> transformations;
  
private:
  DECL_LOGGER("RecoveredInformation");
}; 

template<class graph_pack_t>
void RecoveredInformation<graph_pack_t>::init_raw_results() { 
  for (auto im = graph_pack.multicolors.cbegin_vec_T_consistent_color(); im != graph_pack.multicolors.cend_vec_T_consistent_color(); ++im) {
    graphs.insert(std::make_pair(*im, graph_pack.graph.get_partgraph(0)));
    transformations.insert(std::make_pair(*im, transform_t())); 
  } 
  graphs.insert(std::make_pair(graph_pack.multicolors.get_root_color(), graph_pack.graph.get_partgraph(0)));
  transformations.insert(std::make_pair(graph_pack.multicolors.get_root_color(), transform_t())); 

  for (auto it = graph_pack.history.rbegin(); it != graph_pack.history.rend(); ++it) {
    for (auto im = graph_pack.multicolors.cbegin_vec_T_consistent_color(); im != graph_pack.multicolors.cend_vec_T_consistent_color(); ++im) {
      if (it->get_mcolor().includes(*im)) { 
        it->inverse().apply_single(graphs[*im]);
      }

      if (it->get_mcolor() == *im) {
        transformations[*im].push_back(it->inverse());
      }
    }
  }
}

template<class graph_pack_t>
void RecoveredInformation<graph_pack_t>::init_linearize_results() { 
  INFO("Get history from process graph")
  init_raw_results();

  INFO("Start walk on tree and run algorithm for linearize")
  for (auto const & tree : cfg::get().phylotrees) {
    walk_and_linearize(tree.get_root()); 
  }
  INFO("End walk on tree and run algorithm for linearize")
  
  for (auto const & tree : cfg::get().phylotrees) {
    if (tree.is_phylogenetic_root()) { 
      mcolor_t const & left_color = tree.get_root()->get_left_child()->get_data();
      mcolor_t const & right_color = tree.get_root()->get_right_child()->get_data();  
      for (twobreak_t const & twobreak : transformations[right_color]) { 
        transformations[left_color].push_front(twobreak.inverse());
      }
      transformations.erase(right_color);
      break;
    }
  }
}

template<class graph_pack_t>
void RecoveredInformation<graph_pack_t>::init_target_results() { 
  ;
}

template<class graph_pack_t>
void RecoveredInformation<graph_pack_t>::walk_and_linearize(std::unique_ptr<node_t> const & current) {
  using namespace linearize;
  auto const & left = current->get_left_child(); 
  auto const & right = current->get_right_child(); 

  if (left) { 
    walk_and_linearize(left); 
  }

  if (right) { 
    walk_and_linearize(right); 
  }

  if (left && right && graphs.find(current->get_data()) != graphs.end()) {     
    size_t count_left = count_circular_chromosome(graph_pack, graphs[left->get_data()]); 
    size_t central = count_circular_chromosome(graph_pack, graphs[current->get_data()]); 
    size_t count_right = count_circular_chromosome(graph_pack, graphs[right->get_data()]); 

    if (count_left == 0 && central != 0 && count_right == 0) { 
      transform_t linearize_history = linearize_algo.linearize(graphs[current->get_data()], transformations[left->get_data()], graphs[left->get_data()], central); 
      
      //Apply linearizes twobreaks and modify transformation
      for (twobreak_t const & twobreak : linearize_history) { 
        transformations[current->get_data()].push_back(twobreak);
        twobreak.apply_single(graphs[current->get_data()]); 
        transformations[right->get_data()].push_front(twobreak.inverse());  
      } 
      //transformations[left->get_data()] = new_history.second;

      assert(count_circular_chromosome(graph_pack, graphs[current->get_data()]) == count_circular_chromosome(graph_pack, graphs[left->get_data()]));//Check that all is good
      assert(apply_transformation(graphs[current->get_data()], transformations[left->get_data()]) == graphs[left->get_data()]); //Check left
      assert(apply_transformation(graphs[current->get_data()], transformations[right->get_data()]) == graphs[right->get_data()]); //Check right
      if (current->get_parent() != nullptr) {  
        auto iter = graphs.find(current->get_parent()->get_data());
        if (iter != graphs.end()) {  
          assert(apply_transformation(iter->second, transformations[current->get_data()]) == graphs[current->get_data()]); //Check parrent
        }
      } 
    }
  }
} 

template<class graph_pack_t>
void RecoveredInformation<graph_pack_t>::init_parent_colors(std::unique_ptr<node_t> const & current) {
  auto const & left = current->get_left_child(); 
  auto const & right = current->get_right_child(); 

  if (left) { 
    init_parent_colors(left); 
  }

  if (right) { 
    init_parent_colors(right); 
  }

  if (left && right && graph_pack.multicolors.is_T_consistent_color(current->get_data())) { 
    parent_colors.insert(std::make_pair(left->get_data(), current->get_data()));
    parent_colors.insert(std::make_pair(right->get_data(), current->get_data()));
  } 
}

} 

#endif