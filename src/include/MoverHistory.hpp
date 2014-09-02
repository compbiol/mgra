#ifndef MOVER_HISTORY_HPP
#define MOVER_HISTORY_HPP

template<class graph_t> 
struct MoverHistory {

  typedef typename graph_t::edge_t edge_t;
  typedef typename graph_t::twobreak_t twobreak_t; 
  typedef typename graph_t::partgraph_t partgraph_t; 
  
  MoverHistory(graph_t const & graph, partgraph_t const & b_edges) 
  : m_graph(graph)
  , bad_edges(b_edges) 
  { 
  }

  template<class ClassIterator>
  bool swap_two_break(ClassIterator first, ClassIterator second) const { 
    return first_type_swap_twobreak(first, second);
  }

  template<class ClassIterator>
  bool move_insertion_to_begin(ClassIterator first, ClassIterator second) const {          
    if ((second->get_vertex(0) != Infty && second->get_vertex(1) != Infty 
      && second->get_vertex(0) == this->m_graph.get_obverse_vertex(second->get_vertex(1)) 
      && bad_edges.defined(second->get_vertex(0), second->get_vertex(1))) 
      || (second->get_vertex(2) != Infty && second->get_vertex(3) != Infty 
      && second->get_vertex(2) == this->m_graph.get_obverse_vertex(second->get_vertex(3)) 
      && bad_edges.defined(second->get_vertex(2), second->get_vertex(3)))) {
      return first_type_swap_twobreak(first, second);
    } else { 
      assert(false);
    }
    return false; 
  }

  template<class ClassIterator>
  bool move_insertion_to_end(ClassIterator first, ClassIterator second) const { 
    if (first->get_vertex(0) != Infty && first->get_vertex(1) != Infty 
      && first->get_vertex(0) == this->m_graph.get_obverse_vertex(first->get_vertex(1)) 
      && bad_edges.defined(first->get_vertex(0), first->get_vertex(1))) { 

      return second_type_swap_twobreak(first, second); 

    } else if (first->get_vertex(2) != Infty && first->get_vertex(3) != Infty 
      && first->get_vertex(2) == this->m_graph.get_obverse_vertex(first->get_vertex(3)) 
      && bad_edges.defined(first->get_vertex(2), first->get_vertex(3))) {

      return first_type_swap_twobreak(first, second);

    } else { 
      assert(false);
    }

    return false; 
  }

  template<class ClassIterator>
  bool move_deletion_to_begin(ClassIterator first, ClassIterator second) const { 
    if (second->get_vertex(0) != Infty && second->get_vertex(2) != Infty 
      && second->get_vertex(0) == this->m_graph.get_obverse_vertex(second->get_vertex(2)) 
      && bad_edges.defined(second->get_vertex(0), second->get_vertex(2))) { 

      if ((second->get_arc(0) == std::make_pair(first->get_vertex(0), first->get_vertex(2)))
         || (second->get_arc(0) == std::make_pair(first->get_vertex(1), first->get_vertex(3)))
         || (second->get_arc(1) == std::make_pair(first->get_vertex(0), first->get_vertex(2)))
         || (second->get_arc(1) == std::make_pair(first->get_vertex(1), first->get_vertex(3)))) { 
        return first_type_swap_twobreak(first, second);
      } else { 
        return second_type_swap_twobreak(first, second); 
      }

    } else if (second->get_vertex(1) != Infty && second->get_vertex(3) != Infty 
      && second->get_vertex(1) == this->m_graph.get_obverse_vertex(second->get_vertex(3)) 
      && bad_edges.defined(second->get_vertex(1), second->get_vertex(3))) { 

      if ((second->get_arc(0) == std::make_pair(first->get_vertex(0), first->get_vertex(2)))
         || (second->get_arc(0) == std::make_pair(first->get_vertex(1), first->get_vertex(3)))
         || (second->get_arc(1) == std::make_pair(first->get_vertex(0), first->get_vertex(2)))
         || (second->get_arc(1) == std::make_pair(first->get_vertex(1), first->get_vertex(3)))) { 
        return second_type_swap_twobreak(first, second); 
      } else { 
        return first_type_swap_twobreak(first, second);
      } 

    } else { 
      assert(false);
    }

    return false; 
  }

  template<class ClassIterator>
  bool move_deletion_to_end(ClassIterator first, ClassIterator second) const { 
    if ((first->get_vertex(0) != Infty && first->get_vertex(2) != Infty 
      && first->get_vertex(0) == this->m_graph.get_obverse_vertex(first->get_vertex(2)) 
      && bad_edges.defined(first->get_vertex(0), first->get_vertex(2)))
      || (first->get_vertex(1) != Infty && first->get_vertex(3) != Infty 
      && first->get_vertex(1) == this->m_graph.get_obverse_vertex(first->get_vertex(3)) 
      && bad_edges.defined(first->get_vertex(1), first->get_vertex(3)))) { 
      return first_type_swap_twobreak(first, second);
    } else { 
      assert(false);
    }
    return false; 
  }

private:
  enum swap_t {first_type, second_type};

  template<class ClassIterator>
  bool first_type_swap_twobreak(ClassIterator first, ClassIterator second) const { 
    /*
      If independent replace
      If j_1 and j_2 is weakly dependent then replace
        j_1: (x1,x2) x (y1,y2)
        j_2: (x1,y1) x (x3,y3)      
      into:      
        first_j:  (x1,x2) x (x3,y3)
        second_j: (y3,x2) x (y1,y2)
      If strong dependent - return false
    */
    return swap_iter(first, second, first_type); 
  }

  template<class ClassIterator>
  bool second_type_swap_twobreak(ClassIterator first, ClassIterator second) const { 
    /*
      If independent replace
      If j_1 and j_2 is weakly dependent then replace
        j_1: (x1,x2) x (y1,y2)
        j_2: (x1,y1) x (x3,y3)      
      into:      
        first_j:  (y1,y2) x (y3,x3)
        second_j: (x1,x2) x (x3,y2)
      If strong dependent - return false
    */
    return swap_iter(first, second, second_type);  
  }

  template<class ClassIterator>
  bool swap_iter(ClassIterator first, ClassIterator second, swap_t const & swap_type) const;

private: 
  graph_t const & m_graph;
  partgraph_t const & bad_edges;
};

template<class graph_t> 
template<class ClassIterator>
bool MoverHistory<graph_t>::swap_iter(ClassIterator first, ClassIterator second, swap_t const & swap_type) const { 
  twobreak_t const & j_1 = *first;
  twobreak_t const & j_2 = *second;
  bool result = false; 

  if (j_1.is_dependent(j_2) == 0) { 
    std::iter_swap(first, second);
    result = true; 
  } else if (j_1.is_dependent(j_2) == 1) { 
    edge_t p1;
    edge_t q1;
    edge_t q2;

    auto is_use_edge_lambda = [&] (size_t ind_arc, size_t ind1, size_t ind2, bool is_reverse, bool is_another_reverse) -> bool { 
      bool usage = false; 
      if (j_2.get_arc(ind_arc) == edge_t(j_1.get_vertex(ind1), j_1.get_vertex(ind2))) { 
        usage = true;
  
        if (is_reverse) { 
          p1 = edge_t(j_1.get_arc(0).second, j_1.get_arc(0).first);
          q1 = edge_t(j_1.get_arc(1).second, j_1.get_arc(1).first);
        } else {  
          p1 = j_1.get_arc(0);  
          q1 = j_1.get_arc(1); 
        } 

        if (is_another_reverse) {
          q2 = edge_t(j_2.get_arc(1 - ind_arc).second, j_2.get_arc(1 - ind_arc).first);                    
        } else { 
          q2 = j_2.get_arc(1 - ind_arc);
        }
      }
      return usage; 
    };

    if (swap_type == swap_t::first_type) { 
      bool use_arc = false; 
      for (size_t i = 0; (i < 2) && !use_arc; ++i) {
        use_arc = is_use_edge_lambda(i, 0, 2, false, false) || is_use_edge_lambda(i, 2, 0, false, true)
                    || is_use_edge_lambda(i, 1, 3, true, false) || is_use_edge_lambda(i, 3, 1, true, true);
      }

      if (use_arc) {
        *first = twobreak_t(p1, q2, j_2.get_mcolor());
        *second = twobreak_t(q2.second, p1.second, q1.first, q1.second, j_2.get_mcolor());
        result = true; 
      } 
    } else if (swap_type == swap_t::second_type) { 
      bool use_arc = false; 
      for (size_t i = 0; (i < 2) && !use_arc; ++i) {
        use_arc = is_use_edge_lambda(i, 0, 2, false, true) || is_use_edge_lambda(i, 2, 0, false, false)
                    || is_use_edge_lambda(i, 1, 3, true, true) || is_use_edge_lambda(i, 3, 1, true, false);
      }

      if (use_arc) {
        *first = twobreak_t(q1, q2, j_2.get_mcolor());
        *second = twobreak_t(p1.first, p1.second, q2.second, q1.second, j_2.get_mcolor());
        result = true;
      } 
    } else { 
      assert(false);
    }
  } 

  return result;  
}

/*
for (size_t i = 0; i < 2; ++i) { 
  if (j_2.get_arc(i) == edge_t(j_1.get_vertex(0), j_1.get_vertex(2))) { 
    usearc = true;

    p1 = j_1.get_arc(0); //i 
    q1 = j_1.get_arc(1); //1 - i

    q2 = j_2.get_arc(1 - i);
  } else if (j_2.get_arc(i) == edge_t(j_1.get_vertex(2), j_1.get_vertex(0))) {
    usearc = true;
  
    p1 = j_1.get_arc(0);
    q1 = j_1.get_arc(1);

    q2 = edge_t(j_2.get_arc(1 - i).second, j_2.get_arc(1 - i).first);
  } else if (j_2.get_arc(i) == edge_t(j_1.get_vertex(1), j_1.get_vertex(3))) {
    usearc = true;
  
    p1 = edge_t(j_1.get_arc(0).second, j_1.get_arc(0).first);
    q1 = edge_t(j_1.get_arc(1).second, j_1.get_arc(1).first);
    
    q2 = j_2.get_arc(1 - i);
  } else if (j_2.get_arc(i) == edge_t(j_1.get_vertex(3), j_1.get_vertex(1))) {
    usearc = true;
  
    p1 = edge_t(j_1.get_arc(0).second, j_1.get_arc(0).first);
    q1 = edge_t(j_1.get_arc(1).second, j_1.get_arc(1).first);

    q2 = edge_t(j_2.get_arc(1 - i).second, j_2.get_arc(1 - i).first);
  }
} 
*/
#endif