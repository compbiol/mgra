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

  change_history_t linearize(partgraph_t P, transform_t const & transform, partgraph_t const & Q) const override; 

  void move_insertion_to_begin(transform_t & transformation, citer_transform first, citer_transform second) const {          
  	this->swap_two_twobreaks(transformation, first, second);
  }

  void move_deletion_to_end(transform_t & transformation, citer_transform first, citer_transform second) const { 
  	this->swap_two_twobreaks(transformation, first, second);
  }

  void move_deletion_to_begin(transform_t & transformation, citer_transform first, citer_transform second) const;

  void move_insertion_to_end(transform_t & transformation, citer_transform first, citer_transform second) const; 
}; 


template<class graph_pack_t>
typename IndelLinearize<graph_pack_t>::change_history_t IndelLinearize<graph_pack_t>::linearize(partgraph_t P, transform_t const & transform, partgraph_t const & Q) const { 
	transform_t replace_transformation; transform_t transformation = transform; 

  size_t c_P = count_circular_chromosome(this->graph_pack, P); 
  size_t c_Q = count_circular_chromosome(this->graph_pack, Q);

  while (c_P != c_Q) { 
  	bool changed = false; 
    partgraph_t current = P;
    transform_t removed_chromosome; 

    //Find first twobreak is decrease number of circular chromosomes
    for (auto first_j = transformation.begin(); !changed && first_j != transformation.end();) {     
      first_j->apply_single(current);
      size_t c_PP = count_circular_chromosome(this->graph_pack, current);

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