#ifndef LINEARIZATOR_HPP
#define LINEARIZATOR_HPP

#include "algo/linearization/MoverHistory.hpp"

template<class graph_pack_t>
struct Linearizator {

  using mcolor_t = typename graph_pack_t::mcolor_type;
  using edge_t = typename graph_pack_t::edge_t;
  using twobreak_t = typename graph_pack_t::twobreak_t; 
  using transform_t = typename graph_pack_t::transform_t;
  using partgraph_t = typename graph_pack_t::partgraph_t; 
  
  using history_t = std::tuple<transform_t, transform_t, transform_t>;
  using change_history_t = std::pair<transform_t, transform_t>;

  Linearizator(graph_pack_t const & gp) 
  : graph_pack(gp)
  , mover_history(gp)
  {
  }

  change_history_t linearizate(partgraph_t const & P, transform_t const & history, partgraph_t const & Q) const; 

  size_t count_circular_chromosome(partgraph_t const & local_graph) const;

private:
  transform_t classical_linearization(partgraph_t P, transform_t & transformation, partgraph_t const & Q) const; 

  transform_t deletion_linearization(partgraph_t P, transform_t & transformation, partgraph_t const & Q) const; 

  history_t split_history(transform_t transformation) const;
  
  //GO TO LAMBDA
  bool is_circular_chromosome(partgraph_t const & local_graph, vertex_t const & x, std::unordered_set<vertex_t>& processed) const;

private:
  graph_pack_t const & graph_pack;
  MoverHistory<graph_pack_t> mover_history;
};

template<class graph_pack_t>
typename Linearizator<graph_pack_t>::change_history_t Linearizator<graph_pack_t>::linearizate(partgraph_t const & P, transform_t const & history, partgraph_t const & Q) const { 
  transform_t del_transform; 
  transform_t basic_transform;
  transform_t ins_transform; 

  //std::cerr << "Start split history " << std::endl; 
  std::tie(del_transform, basic_transform, ins_transform) = split_history(history);

  /* calculate P' */
  //std::cerr << "Apply deletion history" << std::endl;
  partgraph_t PP = P; 
  for (twobreak_t const & twobreak : del_transform) { 
    //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
    //   << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
    twobreak.apply_single(PP);
  }

  /* calculate Q' */
  //std::cerr << "Apply basic history" << std::endl;
  partgraph_t QQ = PP;
  for (twobreak_t const & twobreak : basic_transform) { 
    //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
    //   << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
    twobreak.apply_single(QQ);
  }

  /* Test thar split history is good */
  //std::cerr << "Apply insertion history" << std::endl;
  partgraph_t tested = QQ;
  for (twobreak_t const & twobreak : ins_transform) { 
    //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
    //   << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
    twobreak.apply_single(tested);
  }
  assert(tested == Q);
  
  /*Go to algorithm*/
  size_t c_P = count_circular_chromosome(P);
  size_t c_PP = count_circular_chromosome(PP);     
  size_t c_QQ = count_circular_chromosome(QQ); 
  size_t c_Q = count_circular_chromosome(Q);

  assert(c_P > c_Q); 

  //std::cerr << c_P << " " << c_PP << " " << c_QQ << " " << c_Q << std::endl;

  /* Run deletion linearization */
  transform_t del_replace_transform; 
  if (c_P > c_PP) {
    //std::cerr << "Start replace deletion history" << std::endl; 
    del_replace_transform = deletion_linearization(P, del_transform, PP);
    /*for (twobreak_t const & twobreak : del_replace_transform) { 
      std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
         << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
    }*/
    //std::cerr << "End replace deletion history" << std::endl;
  }

  /* Run classical linearization.*/
  transform_t basic_replace_transform; 
  if (c_PP > c_QQ) {
    //std::cerr << "Start replace basic history" << std::endl;
    basic_replace_transform = classical_linearization(PP, basic_transform, QQ);
    //std::cerr << "End replace basic history" << std::endl;
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

  //std::cerr << "End work linearization. We have " << std::endl; 

  return std::make_pair(del_replace_transform, del_transform);
}

template<class graph_pack_t>
typename Linearizator<graph_pack_t>::transform_t Linearizator<graph_pack_t>::classical_linearization(partgraph_t P, transform_t & transformation, partgraph_t const & Q) const {
  transform_t replace_transformation; 

  size_t c_P = count_circular_chromosome(P); 
  size_t c_Q = count_circular_chromosome(Q);

  while (c_P != c_Q) { 
    partgraph_t current = P;
    bool changed = false;
    bool is_all_swap = true;  

    /*Do that first twobreak is decrease number of circular chromosomes*/
    for (auto first_j = transformation.begin(); !changed && first_j != transformation.end();) {     
      first_j->apply_single(current);
      size_t c_PP = count_circular_chromosome(current);

      if (c_P > c_PP) {
        /*Find good twobreak, start to move in begining history*/
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

    /*In head history good twobreak, which c_P > c_P1*/
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

    /*Find first twobreak is decrease number of circular chromosomes*/
    for (auto first_j = transformation.begin(); !changed && first_j != transformation.end();) {     
      //std::cerr << "See on twobreak: " << first_j->get_vertex(0) << " " << first_j->get_vertex(1) << " " << first_j->get_vertex(2) 
      //      << " " << first_j->get_vertex(3) << " " << std::endl;

      first_j->apply_single(current);
      size_t c_PP = count_circular_chromosome(current);

      if (c_P > c_PP) {
        changed = true; 
        
        //std::cerr << "Start to changed this two-break" << std::endl;

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
            bool flag = check_dependent_lambda(*p); 
            if (flag) {
              removed_chromosome.push_front(*p);
              transformation.erase(p++); 
            } 
          }  

          bool flag = check_dependent_lambda(*transformation.begin()); 
          if (flag) {
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
typename Linearizator<graph_pack_t>::history_t Linearizator<graph_pack_t>::split_history(transform_t transformation) const {  
  transform_t deletions; 
  transform_t twobreaks; 
  transform_t insertions;
  
  //std::cerr << "move deletion to the end " << std::endl;
  for (auto twobreak = transformation.begin(); twobreak != transformation.end();) { 
    vertex_t const & p = twobreak->get_vertex(0);
    vertex_t const & q = twobreak->get_vertex(1);
    vertex_t const & x = twobreak->get_vertex(2);
    vertex_t const & y = twobreak->get_vertex(3);

    if ((p != Infty && x != Infty && p == this->graph_pack.graph.get_obverse_vertex(x) && graph_pack.is_prosthetic_chromosome(p, x)) || 
        (q != Infty && y != Infty && q == this->graph_pack.graph.get_obverse_vertex(y) && graph_pack.is_prosthetic_chromosome(q, y))) { 
      bool is_changed = true; 
      auto second_j = twobreak; 

      if (twobreak != transformation.begin()) {  
        for (--twobreak; is_changed && twobreak != transformation.begin(); --twobreak) {
          is_changed = mover_history.move_deletion_to_begin(twobreak, second_j);
          
          if (!is_changed) { 
            //std::cerr << "Strong in deletion" << std::endl;
            transformation.erase(twobreak++);
            transformation.erase(twobreak++);
          } 

          second_j = twobreak; 
        }

        if (is_changed) { 
          is_changed = mover_history.move_deletion_to_begin(transformation.begin(), second_j);
          if (!is_changed) { 
            //std::cerr << "is strong deletion" <<std::endl;
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

  /*Move transformation on end*/
  //std::cerr << "move insertion to the end " << std::endl;
  for (auto twobreak = transformation.begin(); twobreak != transformation.end();) { 
    vertex_t const & p = twobreak->get_vertex(0);
    vertex_t const & q = twobreak->get_vertex(1);
    vertex_t const & x = twobreak->get_vertex(2);
    vertex_t const & y = twobreak->get_vertex(3);

    if ((p != Infty && q != Infty && p == this->graph_pack.graph.get_obverse_vertex(q) && graph_pack.is_prosthetic_chromosome(p, q)) || 
        (x != Infty && y != Infty && x == this->graph_pack.graph.get_obverse_vertex(y) && graph_pack.is_prosthetic_chromosome(x, y))) {     
      bool is_changed = true; 
      auto first_j = twobreak; 

      if (twobreak != (--transformation.end())) {  
        for (++twobreak; is_changed && twobreak != transformation.end(); ++twobreak) {
          is_changed = mover_history.move_insertion_to_end(first_j, twobreak);
          if (!is_changed) { 
            //std::cerr << "Strong in insertion" << std::endl;
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
  
  return history_t(deletions, transformation, insertions);
}

template<class graph_pack_t>
bool Linearizator<graph_pack_t>::is_circular_chromosome(partgraph_t const & local_graph, vertex_t const & x, std::unordered_set<vertex_t>& processed) const {
  bool circular = false;
  bool have_deletion = false;
  vertex_t previous = x;
  vertex_t current = x; 
  processed.insert(x);

  do { 
    previous = current; 
    current = graph_pack.graph.get_obverse_vertex(previous);
    if (processed.count(current) == 0) {
      processed.insert(current);
      if (local_graph.defined(current)) {
        previous = current; 
        current = local_graph[previous];
        if (graph_pack.is_prosthetic_chromosome(previous, current)) {
          have_deletion = true;
        }
        if (processed.count(current) != 0) {
          circular = true;
        } 
        processed.insert(current); 
      }
    } else { 
      circular = true;
    } 
  } while ((current != Infty) && local_graph.defined(current) && !circular);

  if (!circular && local_graph.defined(x)) {  
    for (vertex_t y = local_graph[x]; local_graph.defined(y) && (y != Infty); y = local_graph[y]) {
      processed.insert(y);
      if (y != Infty) {
        y = graph_pack.graph.get_obverse_vertex(y);
        processed.insert(y);
      }
    }
  } 

  if (have_deletion) { 
    return false;
  } 

  return circular;
}

template<class graph_pack_t>
size_t Linearizator<graph_pack_t>::count_circular_chromosome(partgraph_t const & local_graph) const {
  std::unordered_set<vertex_t> processed;
  size_t count_chr = 0;
  
  for (vertex_t const & x : graph_pack.graph) { 
    if (processed.count(x) == 0) { 
      std::unordered_set<vertex_t> chr_set;
      bool circular = is_circular_chromosome(local_graph, x, chr_set); 
      std::copy(chr_set.begin(), chr_set.end(), std::inserter(processed, processed.end()));
      if (circular) {
        ++count_chr; 
      }
    } 
  }
  
  return count_chr;
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