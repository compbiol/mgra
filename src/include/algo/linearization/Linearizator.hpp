#ifndef LINEARIZE_HPP
#define ABS_LINEARIZE_HPP

#include "algo/linearization/MoverHistory.hpp"

namespace algo { 

template<class graph_pack_t>
struct Linearizator {

  using mcolor_t = typename graph_pack_t::mcolor_type;
  using edge_t = typename graph_pack_t::edge_t;
  using twobreak_t = typename graph_pack_t::twobreak_t; 
  using transform_t = typename graph_pack_t::transform_t;
  using partgraph_t = typename graph_pack_t::partgraph_t; 
  
  using genome_t = structure::Genome;
  using chromosome_t = structure::Chromosome; 
  
  using change_history_t = std::pair<transform_t, transform_t>;

  Linearizator(graph_pack_t const & gp) 
  : graph_pack(gp)
  , mover_history(gp)
  {
  }

  /**
   * main procedure for linearization algorithm. See paper for detailed information. 
   * get transformation HISTORY between genome P and genome Q, where c(P) > c(Q). 
   * return pair of transformations T_1, and T_2. 
   * where genome after T_1 have c(P') == c(Q) and length of T_1 is equal c(P) - c(Q). 
   * T_2 is transformation between genome P' -> Q
   */
  change_history_t linearizate(partgraph_t const & P, transform_t const &  history, partgraph_t const & Q) const; 

private:
  using history_t = std::tuple<transform_t, transform_t, transform_t>;

  /**
  * Split input transformation on deletion transformation -> classical transformation -> insertion transformation. 
  * Return tuple with deletions, two-breaks, insertions operations. 
  */
  transform_t get_safely_deletions(transform_t & transformation) const;
  transform_t get_safely_insertions(transform_t & transformation) const;

  transform_t classical_linearization(partgraph_t P, transform_t & transformation, partgraph_t const & Q) const; 
  transform_t deletion_linearization(partgraph_t P, transform_t & transformation, partgraph_t const & Q) const; 

  
private:
  graph_pack_t const & graph_pack;
  MoverHistory<graph_pack_t> mover_history;

private: 
  DECL_LOGGER("LinearizationAlgorithm");
};

template<class graph_pack_t>
typename Linearizator<graph_pack_t>::transform_t Linearizator<graph_pack_t>::get_safely_deletions(transform_t & transformation) const { 
  transform_t deletions; 

  for (auto twobreak = transformation.begin(); twobreak != transformation.end();) { 
    vertex_t const & p = twobreak->get_vertex(0); vertex_t const & q = twobreak->get_vertex(1);
    vertex_t const & x = twobreak->get_vertex(2); vertex_t const & y = twobreak->get_vertex(3);

    if ((p != Infty && x != Infty && p == this->graph_pack.graph.get_obverse_vertex(x) && graph_pack.is_prosthetic_chromosome(p, x)) || 
        (q != Infty && y != Infty && q == this->graph_pack.graph.get_obverse_vertex(y) && graph_pack.is_prosthetic_chromosome(q, y))) { 
      bool is_changed = true; 
      auto second_j = twobreak; 

      if (twobreak != transformation.begin()) {  
        for (--twobreak; is_changed && twobreak != transformation.begin(); --twobreak) {
          is_changed = mover_history.move_deletion_to_begin(twobreak, second_j);
          if (!is_changed) { 
            transformation.erase(twobreak++);
            transformation.erase(twobreak++);
          } 
          second_j = twobreak; 
        }
        if (is_changed) {
          is_changed = mover_history.move_deletion_to_begin(transformation.begin(), second_j);
          if (!is_changed) { 
            transformation.pop_front();
            transformation.pop_front();
          }
        }
      } 
      if (is_changed) { 
        deletions.push_back(transformation.front());    
        transformation.pop_front(); 
      } 
      twobreak = transformation.begin(); 
    } else { 
      ++twobreak;
    }
  }

  return deletions;
}

template<class graph_pack_t>
typename Linearizator<graph_pack_t>::transform_t Linearizator<graph_pack_t>::get_safely_insertions(transform_t & transformation) const { 
  transform_t insertions;

  for (auto twobreak = transformation.begin(); twobreak != transformation.end();) { 
    vertex_t const & p = twobreak->get_vertex(0); vertex_t const & q = twobreak->get_vertex(1);
    vertex_t const & x = twobreak->get_vertex(2); vertex_t const & y = twobreak->get_vertex(3);

    if ((p != Infty && q != Infty && p == this->graph_pack.graph.get_obverse_vertex(q) && graph_pack.is_prosthetic_chromosome(p, q)) || 
        (x != Infty && y != Infty && x == this->graph_pack.graph.get_obverse_vertex(y) && graph_pack.is_prosthetic_chromosome(x, y))) {     
      bool is_changed = true; 
      auto first_j = twobreak; 
      if (twobreak != (--transformation.end())) {  
        for (++twobreak; is_changed && twobreak != transformation.end(); ++twobreak) {
          is_changed = mover_history.move_insertion_to_end(first_j, twobreak);
          if (!is_changed) { 
            transformation.erase(first_j++); 
            transformation.erase(twobreak++); 
          }
          first_j = twobreak; 
        }
      } 
      if (is_changed) { 
        insertions.push_front(transformation.back());
        transformation.pop_back(); 
      } 
      twobreak = transformation.begin(); 
    } else { 
      ++twobreak;
    }
  } 

  return insertions;
}

template<class graph_pack_t>
typename Linearizator<graph_pack_t>::change_history_t Linearizator<graph_pack_t>::linearizate(partgraph_t const & P, transform_t const & history, partgraph_t const & Q) const { 
  transform_t basic_transform = history;
  transform_t del_transform = get_safely_deletions(basic_transform);
  transform_t ins_transform = get_safely_insertions(basic_transform); 

  //calculate P'
  partgraph_t PP = P; 
  for (twobreak_t const & twobreak : del_transform) { 
    //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
    //   << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
    twobreak.apply_single(PP);
  }

  //calculate Q'
  partgraph_t QQ = PP;
  for (twobreak_t const & twobreak : basic_transform) { 
    //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
    //   << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
    twobreak.apply_single(QQ);
  }

  //Test that split history is good 
  partgraph_t tested = QQ;
  for (twobreak_t const & twobreak : ins_transform) { 
    twobreak.apply_single(tested);
  }
  assert(tested == Q);
  
  //Go to algorithm
  size_t c_P = count_circular_chromosome(P);
  size_t c_PP = count_circular_chromosome(PP);     
  size_t c_QQ = count_circular_chromosome(QQ); 
  size_t c_Q = count_circular_chromosome(Q);

  assert(c_P > c_Q); 

  TRACE("Number of circular chromosomes" << c_P << " " << c_PP << " " << c_QQ << " " << c_Q)

  //Run deletion linearization
  transform_t del_replace_transform; 
  if (c_P > c_PP) {
    INFO("Start linearization algortihm on deletion history")
    del_replace_transform = deletion_linearization(P, del_transform, PP);
    INFO("Finish linearization algortihm on deletion history")
    //assert(del_replace_transform.size() == (c_P - c_PP)); 
  }

  //Run classical linearization.
  transform_t basic_replace_transform; 
  if (c_PP > c_QQ) {
    INFO("Start linearization algortihm on classic history")
    basic_replace_transform = classical_linearization(PP, basic_transform, QQ);
    INFO("Finish linearization algortihm on classic history")
    //assert(basic_replace_transform.size() == (c_PP - c_QQ)); 
  } 
  
  while (!basic_replace_transform.empty()) { 
    del_transform.push_back(basic_replace_transform.front());
    basic_replace_transform.pop_front();

    auto first_j = (--del_transform.end()); 
    if (first_j != del_transform.begin()) {
      auto second_j = first_j;
      for (--first_j; first_j != del_transform.begin(); --first_j) {
        mover_history.swap_two_break(first_j, second_j);
        second_j = first_j; 
      } 
      mover_history.swap_two_break(first_j, second_j);
    } 

    del_replace_transform.push_back(del_transform.front()); 
    del_transform.pop_front(); 
  } 

  basic_transform.splice(basic_transform.end(), ins_transform);
  del_transform.splice(del_transform.end(), basic_transform);  

  //assert(del_replace_transform.size() == (c_P - c_Q)); 
  return std::make_pair(del_replace_transform, del_transform);
}

#if 0
first_case_swap_twobreak(first, second) { 
  if (first->is_dependent(*second) == independent_twobreak) { 
    // Theorem 2 
    indepenedent_iter_swap(first, second);
  } else (first->is_dependent(*second) == weakly_dependent_twobreak) {
    // Theorem 3  
    weakly_dependent_iter_swap(first, second);
  } else { 
    // Theorem 3 (Special strong depend) 
  } 
}

second_case_swap_twobreak(first, second, third) { 
  if (first->is_dependent(*second) == independent) { 
    if () { 
      //Theorem 4. if a and b belong to two different chromosomes.
      indepenedent_iter_swap(first, second);
    } else { 
      //Theorem 4. if a and b belong to one chromosome.
      speacial_dependent_iter_swap(second, third);
      speacial_dependent_iter_swap(first, second);
    }
  } else (first->is_dependent(*second) == weakly_dependent) {
    if () { 
      //Theorem 5. if a and b and d belong to two different chromosome and one of circular.
      weakly_dependent_iter_swap(first, second);
    } else if () { 
      //Theorem 5. if a and b and d belong to one chromosome. 
      speacial_dependent_iter_swap(second, third);
      speacial_dependent_iter_swap(first, second);
    } else { 
      //Theorem 5. if a and b and d belong to two different chromosomes and both linear.  
      indepenedent_iter_swap(second, third);
      weakly_dependent_iter_swap(first, second);
    }
  } else { 
  } 
}


//move to begining
void one_step_induction(Iterator start_range, Itertator finish_range, partgraph_t current) { 
  size_t c_q0 = count_circular_chromosome(current);
  auto first_it = start_range; 
  size_t c_q1 = c_q0; 

  //find first operation where c(graph_before) > c(graph_after)
  for (; first_it != finish_range && c_q0 <= c_q1; ++first_it) { 
    c_q0 = c_q1; 
    first_it->apply_single(current);
    c_q1 = count_circular_chromosome(current);
  } 

  if (first_it == start_range) { 
    return;
  }

  size_t c_qq2 = c_q1; //after second_it
  first_it->inverse().apply_single(current);
  auto second_it = first_it; 
  size_t c_qq1 = c_q0; //between first_it and second_it
  (--first_it)->inverse().apply_single(current);
  size_t c_qq0 = count_circular_chromosome(current); //before first_it

  while (first_it != start_range) { 

    if (c_q0 >= c_q1 && c_q1 > c_q2) { 
      //Theorem 2 (independent) and 3 (dependent) by paper. 
      bool is_all_swap = mover_history.swap_two_break(first_it, second_it);    
    } else if (c_q0 < c_q1 && c_q1 > c_q2) { 
      //Theorem 4 (independent) and 5 (dependent) by paper. 
      //here we can start scan to the end again to find c_q3 which c_q2 > c_q3 in one speial case. 
      partgraph_t temp = current; 
      first_it->apply_single(temp); second_it->apply_single(temp);
      auto third_it = second_it; 
      (++third_it)->apply_single(temp);
      one_step_induction(third_it, finish_range, temp); 
      bool is_all_swap = mover_history.swap_two_break(first_it, second_it, third_it);    
    } 

    c_qq2 = c_qq1; //after second_it
    second_it = first_it;
    c_qq1 = c_qq0; //between first_it and second_it
    (--first_it)->inverse().apply_single(current);
    c_qq0 = count_circular_chromosome(current); //before first_it
  } 

  //here once that move to equal start range
  if (c_q0 >= c_q1 && c_q1 > c_q2) { 
    //Theorem 2 (independent) and 3 (dependent) by paper. 
    bool is_all_swap = mover_history.swap_two_break(first_iy, second_it);    
  } else if (c_q0 < c_q1 && c_q1 > c_q2) { 
    //Theorem 4 (independent) and 5 (dependent) by paper. 
    //here we can start scan to the end again to find c_q3 which c_q2 > c_q3 in one speial case. 
    partgraph_t temp = current; 
    first_it->apply_single(temp); second_it->apply_single(temp);
    auto third_it = second_it; 
    (++third_it)->apply_single(temp);
    one_step_induction(third_it, finish_range, temp); 
    bool is_all_swap = mover_history.swap_two_break(first_it, second_it, third_it);    
  } 
}


template<class graph_pack_t>
typename Linearizator<graph_pack_t>::transform_t Linearizator<graph_pack_t>::classical_linearization(partgraph_t P, transform_t & transformation, partgraph_t const & Q) const {
  transform_t replace_transformation; 
  
  for (size_t i = 0; i < (count_circular_chromosome(P) - count_circular_chromosome(Q)); ++i) { 
    one_step_induction(transformation.begin(), transformation.end(), P);
    transformation.begin()->apply_single(P);
    replace_transformation.push_back(*transformation.begin());
    transformation.pop_front();  
  }

  return replace_transformation; 
} 
#endif 

template<class graph_pack_t>
typename Linearizator<graph_pack_t>::transform_t Linearizator<graph_pack_t>::classical_linearization(partgraph_t P, transform_t & transformation, partgraph_t const & Q) const {
  transform_t replace_transformation; 
  size_t c_P = count_circular_chromosome(P); 
  size_t c_Q = count_circular_chromosome(Q);
    
  while (c_P != c_Q) { 
    partgraph_t current = P;
    bool changed = false;
    bool is_all_swap = true;  

    //Do that first twobreak is decrease number of circular chromosomes
    for (auto first_j = transformation.begin(); !changed && first_j != transformation.end();) {     
      first_j->apply_single(current);
      size_t c_PP = count_circular_chromosome(current);

      if (c_P > c_PP) {
        //Find good twobreak, start to move in begining history
        changed = true;

        if (first_j != transformation.begin()) {
          auto second_j = first_j;

          for (--first_j; is_all_swap && first_j != transformation.begin(); --first_j) {
            is_all_swap = mover_history.swap_two_break(first_j, second_j);

            if (!is_all_swap) { 
              transformation.erase(first_j++);
              transformation.erase(first_j++);
            }

            second_j = first_j; 
          } 

          if (is_all_swap) {  
            is_all_swap = mover_history.swap_two_break(first_j, second_j);
            if (!is_all_swap) { 
              transformation.erase(first_j++);
              transformation.erase(first_j++);
            } 
          } 
        } 
      } else { 
        ++first_j; 
      }
    } 

    assert(changed);

    //In head history good twobreak, which c_P > c_P1
    if (is_all_swap) { 
      transformation.begin()->apply_single(P);
      replace_transformation.push_back(*transformation.begin());
      transformation.pop_front(); 
    } 

    c_P = count_circular_chromosome(P);
  } 

  return replace_transformation; 
}

template<class graph_pack_t>
typename Linearizator<graph_pack_t>::transform_t Linearizator<graph_pack_t>::deletion_linearization(partgraph_t P, transform_t & transformation, partgraph_t const & Q) const { 
  transform_t replace_transformation; 
  size_t c_P = count_circular_chromosome(P); 
  size_t c_Q = count_circular_chromosome(Q);

  while (c_P != c_Q) { 
    partgraph_t current = P;
    bool changed = false; 
    transform_t removed_chromosome; 

    //Find first twobreak is decrease number of circular chromosomes
    for (auto first_j = transformation.begin(); !changed && first_j != transformation.end();) {     
      first_j->apply_single(current);
      size_t c_PP = count_circular_chromosome(current);

      if (c_P > c_PP) {
        changed = true; 
        if (first_j == transformation.begin()) { 
          removed_chromosome.push_back(*first_j); 
          transformation.pop_front(); 
        } else {  
          removed_chromosome.push_back(*first_j); 
          transformation.erase(first_j++); 

          auto check_dependent_lambda = [&] (twobreak_t const & tested) -> bool { 
            bool result = false; 
            for (auto br = removed_chromosome.begin(); !result && br != removed_chromosome.end(); ++br) { 
              if (br->is_dependent(tested) != 0) { 
                result = true; 
              }  
            }
            return result;
          }; 

          for (auto p = (--first_j); p != transformation.begin(); --p) { 
            if (check_dependent_lambda(*p)) {
              removed_chromosome.push_front(*p);
              transformation.erase(p++); 
            } 
          }  

          if (check_dependent_lambda(*transformation.begin())) {
            removed_chromosome.push_front(*transformation.begin());
            transformation.pop_front(); 
          }     
        } 
      } else { 
        ++first_j;
      }
    } 

    if (removed_chromosome.size() == 1) { 
      removed_chromosome.front().apply_single(P);
      replace_transformation.push_back(removed_chromosome.front());
    } else if (!removed_chromosome.empty()) { 
      twobreak_t temp = removed_chromosome.front(); 
      twobreak_t fission(temp.get_vertex(0), temp.get_vertex(1), Infty, Infty, temp.get_mcolor()); 
      twobreak_t fusion(temp.get_vertex(0), Infty, temp.get_vertex(1), Infty, temp.get_mcolor());

      transformation.push_front(fusion);
      for (twobreak_t const & twobreak : removed_chromosome) { 
        transformation.push_back(twobreak);
      } 
      fission.apply_single(P);
      replace_transformation.push_back(fission);
    }

    c_P = count_circular_chromosome(P);
  } 

  return replace_transformation;
}

template<class graph_pack_t>
structure::Genome Linearizator<graph_pack_t>::get_genome(std::string const & name, partgraph_t const & local_graph) const { 
  std::string const name_chr("chr");
  std::unordered_set<vertex_t> processed;
  genome_t genome(name); 
  size_t count = 1; 

  for (vertex_t const & x : graph_pack.graph) { 
    if (processed.find(x) == processed.end()) { 
      chromosome_t chromosome = get_chromosome(local_graph, x, processed);
      if (chromosome.size() != 0) {
        genome.insert(name_chr + std::to_string(count++), chromosome);
      } 
    }
  } 

  return genome;
}

template<class graph_pack_t>
size_t Linearizator<graph_pack_t>::count_circular_chromosome(partgraph_t const & local_graph) const {
  std::unordered_set<vertex_t> processed;
  size_t count_circular_chr = 0;
  
  for (vertex_t const & x : graph_pack.graph) { 
    if (processed.count(x) == 0) { 
      chromosome_t chromosome = get_chromosome(local_graph, x, processed); 
      if (chromosome.is_circular()) {
        ++count_circular_chr; 
      }
    } 
  }
  
  return count_circular_chr;
}

template<class graph_pack_t>
structure::Chromosome Linearizator<graph_pack_t>::get_chromosome(partgraph_t const & local_graph, vertex_t const & x, std::unordered_set<vertex_t>& processed) const { 
  std::list<std::pair<vertex_t, int> > path;
  bool circular = false;
  bool have_deletion = false;
  vertex_t previous = x;
  vertex_t current = x; 
  processed.insert(x);

  auto const changer_lambda = [&] (vertex_t const & t, bool flag) -> void {
    if (flag) { 
      (*t.rbegin() == 't')?path.push_back(std::make_pair(t.substr(0, t.size() - 1), -1)):path.push_back(std::make_pair(t.substr(0, t.size() - 1), 1));
    } else { 
      (*t.rbegin() == 't')?path.push_front(std::make_pair(t.substr(0, t.size() - 1), 1)):path.push_front(std::make_pair(t.substr(0, t.size() - 1), -1));  
    }  
  };

  do { 
    previous = current; 
    current = graph_pack.graph.get_obverse_vertex(previous); 
  
    if (processed.count(current) == 0) { 
      processed.insert(current);
      changer_lambda(current, true);

      if (local_graph.defined(current)) { 
        previous = current; 
        current = local_graph.find(previous)->second; 

        if (graph_pack.is_prosthetic_chromosome(previous, current)) {
          have_deletion = true;
        }

        if (processed.count(current) != 0) {
          circular = true;
        } else if (current != Infty) { 
          processed.insert(current);  
        }
      } 
    } else { 
      circular = true;
    } 
  } while ((current != Infty) && local_graph.defined(current) && !circular);

  if (!circular && local_graph.defined(x)) {  
    for (vertex_t y = local_graph.find(x)->second; local_graph.defined(y) && (y != Infty); y = local_graph.find(y)->second) {
      processed.insert(y);
      if (y != Infty) {
        y = graph_pack.graph.get_obverse_vertex(y);
        processed.insert(y);
        changer_lambda(y, false);
      }
    } 
  } 
    
  if (have_deletion) {
    return chromosome_t();
  }

  return chromosome_t(path, circular);
}

} 
/*
for (auto q = transformation.begin(); q != transformation.end(); ++q) { 
    q->apply_single(P);
    size_t c_PP = count_circular_chromosome(P); 

    if (c_P > c_PP) {
      transform_t removed_chromosome; 
      removed_chromosome.push_back(q); 

      for (auto p = other_transformation.rbegin(); p != other_transformation.rend();) { 

        bool flag = false; 
        for (auto br = removed_chromosome.begin(); !flag && br != removed_chromosome.end(); ++br) { 
          if (p->is_dependent(*br) != 0) { 
            flag = true; 
            removed_chromosome.push_front(*p);
          }  
        }

        if (flag) {
          other_transformation.erase(std::next(p).base()); 
        } else { 
          ++p;
        }

      }       

      for (twobreak_t const & twobreak : removed_chromosome) {  
        replace_transformation.push_back(twobreak);
      } 
    } else { 
      other_transformation.push_back(*q);
    }
  }

  return std::make_pair(replace_transformation, other_transformation); 
*/

      /*std::cerr << "Swap pair" << second_j->is_dependent(*twobreak) << std::endl;
          std::cerr << "First: " << twobreak->get_vertex(0) << " " << twobreak->get_vertex(1) << " " << twobreak->get_vertex(2) 
            << " " << twobreak->get_vertex(3) << " " << std::endl;
          std::cerr << "Second: " << second_j->get_vertex(0) << " " << second_j->get_vertex(1) << " " << second_j->get_vertex(2) 
            << " " << second_j->get_vertex(3) << " " << std::endl;
          */
      

          /*std::cerr << "After swap " << std::endl;
          std::cerr << "First: " << twobreak->get_vertex(0) << " " << twobreak->get_vertex(1) << " " << twobreak->get_vertex(2) 
            << " " << twobreak->get_vertex(3) << " " << std::endl;
          std::cerr << "Second: " << second_j->get_vertex(0) << " " << second_j->get_vertex(1) << " " << second_j->get_vertex(2) 
            << " " << second_j->get_vertex(3) << " " << std::endl;
          */

#endif