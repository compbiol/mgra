#ifndef FISSION_DEL_LINEARIZE_HPP
#define FISSION_DEL_LINEARIZE_HPP

namespace algo { 

namespace linearize { 

template<class graph_pack_t>
struct IndelLinearize : public algo::linearize::AbsLinearize<graph_pack_t> {
  using mcolor_t = typename graph_pack_t::mcolor_type;
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

  change_history_t linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q) const override; 

  void move_insertion_to_begin(transform_t & transformation, citer_transform first, citer_transform second) const {          
  	this->swap_two_twobreaks(transformation, first, second);
  }

  void move_deletion_to_end(transform_t & transformation, citer_transform first, citer_transform second) const { 
  	this->swap_two_twobreaks(transformation, first, second);
  }

  void move_deletion_to_begin(transform_t & transformation, citer_transform first, citer_transform second) const;

  void move_insertion_to_end(transform_t & transformation, citer_transform first, citer_transform second) const; 

private: 
  std::pair<citer_transform, citer_transform> find_range(citer_transform start, partgraph_t current, citer_transform finish) const {
    size_t c_p = count_circular_chromosome(this->graph_pack, current); size_t c_q = c_p;

    std::pair<citer_transform, citer_transform> range(start, start);
    for (; range.second != finish && (c_p <= c_q); ++range.second) { 
      range.second->apply_single(current);
      c_q = count_circular_chromosome(this->graph_pack, current);
    } 

    return range;
  } 

}; 


template<class graph_pack_t>
typename ClassicalLinearize<graph_pack_t>::citer_transform ClassicalLinearize<graph_pack_t>::one_step_induction(citer_transform start, partgraph_t current, citer_transform finish) const {  
  auto range = find_range(start, current, finish);
  auto last_twobreak = range.second;

  if (finish != range.second) { 

    if (range.second != start) { 
      for (auto iter = (--range.second); iter != start; --iter) { 
        if (iter->is_dependent(*last_twobeak) != twobreak_t::independent) {
         last_twobreak = iter;
        } 
      }  

      if (start->is_dependent(*last_twobeak) != twobreak_t::independent) {
        last_twobreak = start;
      } 
    }

    twobreak_t fission(last_twobreak->get_vertex(0), last_twobreak->get_vertex(1), Infty, Infty, last_twobreak->get_mcolor()); 
    twobreak_t fusion(last_twobreak->get_vertex(0), Infty, last_twobreak->get_vertex(1), Infty, last_twobreak->get_mcolor());
  } 

} 

template<class graph_pack_t>
typename IndelLinearize<graph_pack_t>::change_history_t IndelLinearize<graph_pack_t>::linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q) const {
  auto start = transform.begin(); auto finish = transform.end();

  size_t diff_chromosomes = count_circular_chromosome(this->graph_pack, P) - count_circular_chromosome(this->graph_pack, Q);
  for (size_t i = 0; i < diff_chromosomes; ++i) { 
    std::cerr << "Start step induction " << count_circular_chromosome(this->graph_pack, P) << std::endl;     
    auto spliters = one_step_induction(start, P, finish);
    transform.push_front(spliters.second); transform.push_front(spliters.first);
    start = transform.begin(); finish = transform.end();
    std::cerr << "Finish step induction " << count_circular_chromosome(this->graph_pack, P) << std::endl; 
  }

} 

template<class graph_pack_t>
typename IndelLinearize<graph_pack_t>::change_history_t IndelLinearize<graph_pack_t>::linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q) const { 
	transform_t replace_transformation; transform_t transformation = transform; 
  size_t c_P = count_circular_chromosome(this->graph_pack, P); 
  size_t c_Q = count_circular_chromosome(this->graph_pack, Q);
  twobreak_t last_twobreak; 

  while (c_P != c_Q) { 
    partgraph_t current = P;

    //Find first twobreak is decrease number of circular chromosomes
    auto first_j = transformation.begin()
    size_t c_PP = c_P; 
    for (; first_j != transformation.end() && (c_P <= c_PP); ++first_j) {     
      first_j->apply_single(current);
      size_t c_PP = count_circular_chromosome(this->graph_pack, current);
    } 

    if (transformation.end() != first_j) { 
      last_twobeak = *first_j; 

      if (first_j != transformation.begin()) { 
        for (auto p = (--first_j); p != transformation.begin(); --p) { 
          if (p->is_dependent(last_twobeak) != twobreak_t::independent) {
           last_twobreak = *p;
          } 
        }  

        if (transformation.begin()->is_dependent(last_twobeak) != twobreak_t::independent) {
          last_twobreak = *transformation.begin();
        } 
      } 
    } 

    if (!removed_chromosome.empty()) { 
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

    c_P = count_circular_chromosome(this->graph_pack, P);
  } 

  return std::make_pair(replace_transformation, transformation);
}


template<class graph_pack_t>
void IndelLinearize<graph_pack_t>::move_deletion_to_begin(transform_t & transformation, citer_transform first, citer_transform second) const { 
  if (second->get_vertex(0) != Infty && second->get_vertex(2) != Infty 
    && second->get_vertex(0) == this->graph_pack.graph.get_obverse_vertex(second->get_vertex(2))) { 
    //&& this->graph_pack.is_prosthetic_chromosome(second->get_vertex(0), second->get_vertex(2))) { 

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
    //&& this->graph_pack.is_prosthetic_chromosome(second->get_vertex(1), second->get_vertex(3))) { 

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
    //&& this->graph_pack.is_prosthetic_chromosome(first->get_vertex(0), first->get_vertex(1))) { 
  	this->swap_two_twobreaks(transformation, first, second, AbsLinearize<graph_pack_t>::second_type);
  } else if (first->get_vertex(2) != Infty && first->get_vertex(3) != Infty 
    && first->get_vertex(2) == this->graph_pack.graph.get_obverse_vertex(first->get_vertex(3))) { 
    //&& this->graph_pack.is_prosthetic_chromosome(first->get_vertex(2), first->get_vertex(3))) {
    this->swap_two_twobreaks(transformation, first, second);
  } 
}

}

} 

#endif 