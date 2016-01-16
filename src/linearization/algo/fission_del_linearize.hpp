#ifndef FISSION_DEL_LINEARIZE_HPP
#define FISSION_DEL_LINEARIZE_HPP

namespace algo { 

namespace linearize { 

template<class graph_pack_t>
struct IndelLinearize : public algo::linearize::AbsLinearize<graph_pack_t> {
  using mcolor_t = typename graph_pack_t::mcolor_t;
  using edge_t = typename graph_pack_t::edge_t;
  using twobreak_t = typename graph_pack_t::twobreak_t; 
  using transform_t = typename graph_pack_t::transform_t;
  using partgraph_t = typename graph_pack_t::partgraph_t; 
  using citer_transform = typename transform_t::iterator; 
  using genome_t = structure::Genome;
  using chromosome_t = structure::Chromosome;   
  using change_history_t = std::pair<transform_t, transform_t>;

  explicit IndelLinearize(graph_pack_t const & gp) 
  : AbsLinearize<graph_pack_t>(gp)
  {
  }

  transform_t linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q, size_t diff_chromosomes) const override; 

  void move_insertion_to_begin(transform_t & transformation, citer_transform first, citer_transform second) const {          
  	this->swap_two_twobreaks(transformation, first, second);
  }
  
  void move_deletion_to_end(transform_t & transformation, citer_transform first, citer_transform second) const { 
  	this->swap_two_twobreaks(transformation, first, second);
  }

  void move_deletion_to_begin(transform_t & transformation, citer_transform first, citer_transform second) const;
  void move_insertion_to_end(transform_t & transformation, citer_transform first, citer_transform second) const; 

private: 
  /**
   * Find big ranges [....>]
   */
  std::pair<citer_transform, citer_transform> find_range(citer_transform start, partgraph_t current, citer_transform finish) const {
    size_t c_p = count_circular_chromosome(this->graph_pack, current); size_t c_q = c_p;

    std::pair<citer_transform, citer_transform> range(start, start);
    for (; range.second != finish && (c_p <= c_q); ++range.second) { 
      range.second->apply_single(current);
      c_q = count_circular_chromosome(this->graph_pack, current);
    } 

    return range;
  } 

  typename std::pair<twobreak_t, twobreak_t> one_step_induction(citer_transform start, partgraph_t current, citer_transform finish) const;
}; 

template<class graph_pack_t>
void IndelLinearize<graph_pack_t>::move_deletion_to_begin(transform_t & transformation, citer_transform first, citer_transform second) const { 
  if (second->get_vertex(0) != Infty && second->get_vertex(2) != Infty 
    && second->get_vertex(0) == this->graph_pack.graph.get_obverse_vertex(second->get_vertex(2))) { 
    
    if ((second->get_arc(0) == std::make_pair(first->get_vertex(0), first->get_vertex(2)))
       || (second->get_arc(0) == std::make_pair(first->get_vertex(1), first->get_vertex(3)))
       || (second->get_arc(1) == std::make_pair(first->get_vertex(0), first->get_vertex(2)))
       || (second->get_arc(1) == std::make_pair(first->get_vertex(1), first->get_vertex(3)))) { 
      this->swap_two_twobreaks(transformation, first, second);
    } else { 
      this->swap_two_twobreaks(transformation, first, second, AbsLinearize<graph_pack_t>::second_type);
    }

  } else if (second->get_vertex(1) != Infty && second->get_vertex(3) != Infty 
    && second->get_vertex(1) == this->graph_pack.graph.get_obverse_vertex(second->get_vertex(3))) { 
    
    if ((second->get_arc(0) == std::make_pair(first->get_vertex(0), first->get_vertex(2)))
       || (second->get_arc(0) == std::make_pair(first->get_vertex(1), first->get_vertex(3)))
       || (second->get_arc(1) == std::make_pair(first->get_vertex(0), first->get_vertex(2)))
       || (second->get_arc(1) == std::make_pair(first->get_vertex(1), first->get_vertex(3)))) { 
      this->swap_two_twobreaks(transformation, first, second, AbsLinearize<graph_pack_t>::second_type);
    } else { 
      this->swap_two_twobreaks(transformation, first, second);
    } 

  } 

}

template<class graph_pack_t>
void IndelLinearize<graph_pack_t>::move_insertion_to_end(transform_t & transformation, citer_transform first, citer_transform second) const { 
  if (first->get_vertex(0) != Infty && first->get_vertex(1) != Infty 
    && first->get_vertex(0) == this->graph_pack.graph.get_obverse_vertex(first->get_vertex(1))) { 
    this->swap_two_twobreaks(transformation, first, second, AbsLinearize<graph_pack_t>::second_type);
  } else if (first->get_vertex(2) != Infty && first->get_vertex(3) != Infty 
    && first->get_vertex(2) == this->graph_pack.graph.get_obverse_vertex(first->get_vertex(3))) { 
    this->swap_two_twobreaks(transformation, first, second);
  } 
}

template<class graph_pack_t>
typename std::pair<typename graph_pack_t::twobreak_t, typename graph_pack_t::twobreak_t> IndelLinearize<graph_pack_t>::one_step_induction(citer_transform start, partgraph_t current, citer_transform finish) const {  
  auto range = find_range(start, current, finish);
  auto last_twobreak = range.second;
  --last_twobreak;
  twobreak_t fission(last_twobreak->get_vertex(0), last_twobreak->get_vertex(1), Infty, Infty, last_twobreak->get_mcolor()); 
  twobreak_t fusion(last_twobreak->get_vertex(0), Infty, last_twobreak->get_vertex(1), Infty, last_twobreak->get_mcolor()); 
  return std::make_pair(fission, fusion);
} 


template<class graph_pack_t>
typename graph_pack_t::transform_t IndelLinearize<graph_pack_t>::linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q, size_t diff_chromosomes) const {
  transform_t linearize_transformation; 
  diff_chromosomes = std::min(count_circular_chromosome(this->graph_pack, P) - count_circular_chromosome(this->graph_pack, Q), diff_chromosomes);

  for (size_t i = 0; i < diff_chromosomes; ++i) { 
    auto spliters = one_step_induction(transform.begin(), P, transform.end());
    linearize_transformation.push_back(spliters.first); 
    spliters.first.apply_single(P);
    transform.push_front(spliters.second);
  }

  return linearize_transformation;
} 

}

} 

#endif 