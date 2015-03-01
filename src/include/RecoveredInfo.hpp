#ifndef RECOVERED_INFO_HPP
#define RECOVERED_INFO_HPP

#include "Linearizator.hpp"

template<class graph_pack_t>
struct RecoveredInfo {
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using twobreak_t = typename graph_pack_t::twobreak_t; 
  using transform_t = typename graph_pack_t::transform_t;
  using partgraph_t = typename graph_pack_t::partgraph_t; 

  using genome_t = structure::Genome;
  using chromosome_t = structure::Chromosome; 
  using history_t = std::tuple<transform_t, transform_t, transform_t>;

  RecoveredInfo(graph_pack_t const & graph_pack);

  std::vector<genome_t> get_genomes() const { 
    std::vector<genome_t> genomes;
     for (auto const & local_graph : graphs) {
      genomes.push_back(get_genome(cfg::get().mcolor_to_name(local_graph.first), local_graph.second)); 
    }
    return genomes;   
  }
  
  std::map<std::pair<mcolor_t, mcolor_t>, transform_t> get_history() const { 
    std::map<std::pair<mcolor_t, mcolor_t>, transform_t> result_transformations;
    for (auto const & elem : transformations) { 
      assert(parent_colors.find(elem.first) != parent_colors.cend());
      result_transformations.insert(std::make_pair(std::make_pair(elem.first, parent_colors.find(elem.first)->second), elem.second));  
    }
    return result_transformations;
  }

private: 
  using tree_t = typename structure::BinaryTree<mcolor_t>; 
  using node_t = typename tree_t::Node; 
  void get_ugly_history();
  void walk_and_linearizeate(std::unique_ptr<node_t> const & current);
  
  genome_t get_genome(std::string const & name, partgraph_t const & recovered_graph) const;
  chromosome_t get_chromosome(partgraph_t const & recovered_graph, vertex_t const & x, std::unordered_set<vertex_t>& chromosome_set) const;

private: 
  graph_pack_t const & graph_pack;
  Linearizator<graph_pack_t> linearizator;

  std::map<mcolor_t, partgraph_t> graphs;
  std::map<mcolor_t, mcolor_t> parent_colors; 
  std::map<mcolor_t, transform_t> transformations;

  
  partgraph_t bad_edges; //FIXME REMOVE
private:
  DECL_LOGGER("RecoveredInfo");
}; 

template<class graph_pack_t>
void RecoveredInfo<graph_pack_t>::get_ugly_history() { 
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
RecoveredInfo<graph_pack_t>::RecoveredInfo(graph_pack_t const & gp) 
: graph_pack(gp)
, linearizator(gp)
, bad_edges(gp.get_bad_edges())
{ 
  assert(cfg::get().how_build == default_algo);

  INFO("Get history from process graph")
  get_ugly_history();
  
  INFO("Start walk on tree and run algorithm for linearizeate")
  for (auto const & tree : cfg::get().phylotrees) {
    walk_and_linearizeate(tree.get_root()); 
  }
  INFO("End walk on tree and run algorithm for linearizeate")
  
  for (auto const & tree : cfg::get().phylotrees) {
    if (tree.is_phylogenetic_root()) { 
      mcolor_t const & left_color = tree.get_root()->get_left_child()->get_data();
      mcolor_t const & right_color = tree.get_root()->get_right_child()->get_data();  
      parent_colors.insert(std::make_pair(left_color, right_color));
      for (twobreak_t const & twobreak : transformations[right_color]) { 
        transformations[left_color].push_front(twobreak.inverse());
      }
      transformations.erase(right_color);
      break;
    }
  }
}

template<class graph_pack_t>
void RecoveredInfo<graph_pack_t>::walk_and_linearizeate(std::unique_ptr<node_t> const & current) {
  auto const & left = current->get_left_child(); 
  auto const & right = current->get_right_child(); 

  if (left) { 
    walk_and_linearizeate(left); 
  }

  if (right) { 
    walk_and_linearizeate(right); 
  }

  if (left && right && graphs.find(current->get_data()) != graphs.end()) { 
    parent_colors.insert(std::make_pair(left->get_data(), current->get_data()));
    parent_colors.insert(std::make_pair(right->get_data(), current->get_data()));
    
    size_t count_left = linearizator.count_circular_chromosome(graphs[left->get_data()]); 
    size_t central = linearizator.count_circular_chromosome(graphs[current->get_data()]); 
    size_t count_right = linearizator.count_circular_chromosome(graphs[right->get_data()]); 

    if (count_left == 0 && central != 0 && count_right == 0) { 
      std::pair<transform_t, transform_t> new_history = linearizator.linearizate(graphs[current->get_data()], transformations[left->get_data()], graphs[left->get_data()]); 

      //Apply linearization twobreaks and modify transformation
      for (twobreak_t const & twobreak : new_history.first) { 
        transformations[current->get_data()].push_back(twobreak);
        twobreak.apply_single(graphs[current->get_data()]); 
        transformations[right->get_data()].push_front(twobreak.inverse());  
      } 
      transformations[left->get_data()] = new_history.second;

      //Check that all is good
      //Check left
      assert(linearizator.count_circular_chromosome(graphs[current->get_data()]) == 0);
      partgraph_t traverse_graph = graphs[current->get_data()]; 
      for (twobreak_t const & twobreak : transformations[left->get_data()]) { 
        twobreak.apply_single(traverse_graph);
      } 
      assert(traverse_graph == graphs[left->get_data()]);

      //Check right
      traverse_graph = graphs[current->get_data()]; 
      for (twobreak_t const & twobreak : transformations[right->get_data()]) { 
        twobreak.apply_single(traverse_graph);
      } 
      assert(traverse_graph == graphs[right->get_data()]);

      //Check parent
      if (current->get_parent() != nullptr) {  
        auto iter = graphs.find(current->get_parent()->get_data());
        if (iter != graphs.end()) {  
          traverse_graph = iter->second; 
          for (twobreak_t const & twobreak : transformations[current->get_data()]) { 
            twobreak.apply_single(traverse_graph);
          } 
          assert(traverse_graph == graphs[current->get_data()]);
        }
      } 
    }
  }
} 

template<class graph_pack_t>
structure::Genome RecoveredInfo<graph_pack_t>::get_genome(std::string const & name, partgraph_t const & recovered_graph) const { 
  genome_t genome(name); 
  
  std::unordered_set<vertex_t> processed;
  std::string name_chr("chr");
  size_t count = 1; 

  for(vertex_t const & x : graph_pack.graph) { 
    if (processed.find(x) == processed.end()) { 
      std::unordered_set<vertex_t> chromosome_set;
      chromosome_t chromosome = get_chromosome(recovered_graph, x, chromosome_set);
      if (chromosome.size() != 0) {
        genome.insert(name_chr + std::to_string(count++), chromosome);
      } 
      std::copy(chromosome_set.cbegin(), chromosome_set.cend(), std::inserter(processed, processed.end()));
    }
  } 

  return genome;
}

template<class graph_pack_t>
structure::Chromosome RecoveredInfo<graph_pack_t>::get_chromosome(partgraph_t const & recovered_graph, vertex_t const & x, std::unordered_set<vertex_t>& chromosome_set) const {
  std::list<std::pair<vertex_t, int> > path;
  bool circular = false;
  bool have_deletion = false;

  auto const changer_lambda = [&] (vertex_t const & t, bool flag) -> void {
    if (flag) { 
      (*t.rbegin() == 't')?path.push_back(std::make_pair(t.substr(0, t.size() - 1), -1)):path.push_back(std::make_pair(t.substr(0, t.size() - 1), 1));
    } else { 
      (*t.rbegin() == 't')?path.push_front(std::make_pair(t.substr(0, t.size() - 1), 1)):path.push_front(std::make_pair(t.substr(0, t.size() - 1), -1));  
    }  
  };

  vertex_t current = graph_pack.graph.get_obverse_vertex(x); 
  vertex_t previous = x;
  chromosome_set.insert(x);
  
  while ((current != Infty) && !circular) { 
    if (chromosome_set.count(current) != 0) {
      circular = true;
    } else { 
      chromosome_set.insert(current);
      changer_lambda(current, true);

      if (!recovered_graph.defined(current)) { 
        break; // linear
      } 

      previous = current;
      current = recovered_graph.find(previous)->second; 

      if (bad_edges.defined(previous, current)) {
        have_deletion = true;
      }

      if (chromosome_set.count(current) != 0) {
        circular = true;
      } else { 
        chromosome_set.insert(current);
        if (current != Infty) { 
          previous = current;
          current = graph_pack.graph.get_obverse_vertex(previous); 
        } 
      } 
    }
  }

  if (!circular && recovered_graph.defined(x)) {  
    vertex_t y = x;
    while (recovered_graph.defined(y) && (y != Infty)) {
      y = recovered_graph.find(y)->second;
      chromosome_set.insert(y);

      if (y != Infty) {
        y = graph_pack.graph.get_obverse_vertex(y);
        chromosome_set.insert(y);
        changer_lambda(y, false);
      }
    }
  } 
    
  if (have_deletion) {
    return chromosome_t();
  }

  return chromosome_t(path, circular);
}

#if 0 
template<class linearizator_t>
void walk_and_linearizeate() const {
  root->walk_and_linearizeate(linearizator, parent_colors, graphs, transformations);

  mcolor_t left = parent_colors[root->get_left_child()->get_data()];
  mcolor_t right = parent_colors[root->get_right_child()->get_data()];

  if (transformations[left].size() < transformations[right].size()) { 
    parent_colors.erase(left); 
    for (auto const & twobreak : transformations[left]) { 
      transformations[right].push_front(twobreak.inverse());
    }
    transformations.erase(left);
  } else { 
    parent_colors.erase(right); 
    for (auto const & twobreak : transformations[right]) { 
      transformations[left].push_front(twobreak.inverse());
    }
    transformations.erase(right);
  }
}

//ERROR, FIXME, think about edge from root in this function
template<class graph_t>
void RecoveredInfo<graph_t>::walk_and_linearizeate() const { 
  if (left_child) { 
    left_child->walk_and_linearizeate(); 
  }

  if (right_child) { 
    right_child->walk_and_linearizeate(); 
  }

  if (left_child && right_child && graphs.find(this->data) != graphs.end()) { 
    size_t count_left = linearizator.count_circular_chromosome(graphs[left_child->data]); 
    size_t central = linearizator.count_circular_chromosome(graphs[this->data]); 
    size_t count_right = linearizator.count_circular_chromosome(graphs[right_child->data]); 
    //std::cerr << "Left have " << count_left << std::endl 
    //      << " Central have " << central << std::endl << "Right have " << count_right << std::endl;

    if (count_left == 0 && count_right == 0 && central != 0) { 
      typedef typename linearizator_t::twobreak_t twobreak_t;
      typedef typename linearizator_t::transform_t transform_t;
      typedef typename linearizator_t::partgraph_t partgraph_t;

      //std::cerr << "Start linearizator " << genome_match::mcolor_to_name(this->data) 
      //  << " -> " << genome_match::mcolor_to_name(left_child->data) << std::endl;
      std::pair<transform_t, transform_t> new_history = linearizator.linearizate(graphs[this->data], transformations[left_child->data], graphs[left_child->data]); 

      /*Apply linearization twobreaks*/      
      for (twobreak_t const & twobreak : new_history.first) {
        //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
        //<< " " << twobreak.get_vertex(3) << p" " << std::endl;    
        twobreak.apply_single(graphs[this->data]); 
        //std::cerr << linearizator.count_circular_chromosome(graphs[this->data]) << std::endl;
      }
      
      /*Check that all is good*/
      //std::cerr << "Result genome have " << linearizator.count_circular_chromosome(graphs[this->data]) << std::endl;
      assert(linearizator.count_circular_chromosome(graphs[this->data]) == 0);

      //std::cerr << "Start to check linearizator " << std::endl;
      partgraph_t current = graphs[this->data]; 
      for (twobreak_t const & twobreak : new_history.second) { 
        //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
        //<< " " << twobreak.get_vertex(3) << " " << std::endl;    
        
        twobreak.apply_single(current);
      } 
      //std::cerr << "Check that we get good history " << std::endl;
      assert(current == graphs[left_child->data]);

      /*Modify transformation*/
      //std::cerr << "modify transformation" << std::endl;
      for (twobreak_t const & twobreak : new_history.first) { 
        transformations[this->data].push_back(twobreak);
        transformations[right_child->data].push_front(twobreak.inverse());  
      } 
      transformations[left_child->data] = new_history.second;

      //std::cerr << "apply on right child" << std::endl;
      current = graphs[this->data]; 
      for (twobreak_t const & twobreak : transformations[right_child->data]) { 
        twobreak.apply_single(current);
      } 
      assert(current == graphs[right_child->data]);

      if (this->parent->parent != nullptr) {  
        current = graphs[this->parent->data]; 
        for (twobreak_t const & twobreak : transformations[this->data]) { 
          twobreak.apply_single(current);
        } 
        assert(current == graphs[this->data]);
      } 
    }
  } 
}
//std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << p" " << std::endl;    
#endif

#endif