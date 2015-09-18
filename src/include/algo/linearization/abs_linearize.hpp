#ifndef ABS_LINEARIZE_HPP
#define ABS_LINEARIZE_HPP

namespace algo { 

namespace linearize { 

template<class graph_pack_t>
struct AbsLinearize {
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using edge_t = typename graph_pack_t::edge_t;
  using twobreak_t = typename graph_pack_t::twobreak_t; 
  using transform_t = typename graph_pack_t::transform_t;
  using partgraph_t = typename graph_pack_t::partgraph_t; 
  using citer_transform = typename transform_t::iterator; 
  using change_history_t = std::pair<transform_t, transform_t>;
  
  enum swap_t {first_type, second_type};

  explicit AbsLinearize(graph_pack_t const & gp) 
  : graph_pack(gp)
  {
  }  

  virtual transform_t linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q, size_t diff_chromosomes) const = 0; 

  virtual ~AbsLinearize() { 
  } 

protected: 
  citer_transform swap_two_twobreaks(citer_transform first, citer_transform second, citer_transform finish, swap_t const & swap_type = first_type) const { 
    if (first->is_dependent(*second) == twobreak_t::independent) { 
      swap_independent_two_twobreaks(first, second); 
    } else if (first->is_dependent(*second) == twobreak_t::weakly_dependent) { 
      swap_transposition_two_twobreaks(first, second, swap_type);
    } else if (first->is_dependent(*second) == twobreak_t::strong_dependent) { 
      finish = swap_strong_dependent_two_twobreaks(first, second, finish);
    }
    return finish;
  }


  void swap_two_twobreaks(transform_t & transformation, citer_transform first, citer_transform second, swap_t const & swap_type = first_type) const { 
    if (first->is_dependent(*second) == twobreak_t::independent) { 
      swap_independent_two_twobreaks(first, second); 
    } else if (first->is_dependent(*second) == twobreak_t::weakly_dependent) { 
      swap_transposition_two_twobreaks(first, second, swap_type);
    } else if (first->is_dependent(*second) == twobreak_t::strong_dependent) { 
      swap_strong_dependent_two_twobreaks(transformation, first, second);
    }
  }

  void swap_independent_two_twobreaks(citer_transform first, citer_transform second) const { 
    std::iter_swap(first, second); 
  }

  /**
   * If j_1 and j_2 is weakly dependent then replace
   *   j_1: (x1,x2) x (y1,y2)
   *   j_2: (x1,y1) x (x3,y3)      
   * into:      
   *   first_j:  (x1,x2) x (x3,y3)
   *   second_j: (y3,x2) x (y1,y2)
   */
  void swap_first_type_weakly_two_twobreaks(citer_transform first, citer_transform second) const { 
    swap_transposition_two_twobreaks(first, second, first_type); 
  }

  /**
   * If j_1 and j_2 is weakly dependent then replace
   *   j_1: (x1,x2) x (y1,y2)
   *   j_2: (x1,y1) x (x3,y3)      
   * into:      
   *  first_j:  (y1,y2) x (y3,x3)
   *  second_j: (x1,x2) x (x3,y2)
   */
  void swap_second_type_weakly_two_twobreaks(citer_transform first, citer_transform second) const { 
    swap_transposition_two_twobreaks(first, second, second_type);  
  }

  citer_transform swap_strong_dependent_two_twobreaks(citer_transform first, citer_transform second, citer_transform finish) const { 
    for (++second; second != finish; ++first, ++second) { 
      *first =  *second;
    } 
    --finish; 
    return (--finish);
  }

  void swap_strong_dependent_two_twobreaks(transform_t & transformation, citer_transform first, citer_transform second) const { 
    std::cerr << "ATTENTION: STRONG DEPENDENT" << std::endl; 
    transformation.erase(first++); 
    transformation.erase(first++);
    second = first; 
  }

private: 
  template<class ClassIterator>
  void swap_transposition_two_twobreaks(ClassIterator first, ClassIterator second, swap_t const & swap_type) const;

protected:
  graph_pack_t const & graph_pack;
};

template<class graph_pack_t> 
template<class ClassIterator>
void AbsLinearize<graph_pack_t>::swap_transposition_two_twobreaks(ClassIterator first, ClassIterator second, swap_t const & swap_type) const { 
  twobreak_t const & j_1 = *first;
  twobreak_t const & j_2 = *second;
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
    } 
  } 
} 

} 

} 

#endif 