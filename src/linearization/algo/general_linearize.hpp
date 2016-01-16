#ifndef GENERAL_LINEARIZE_HPP
#define GENERAL_LINEARIZE_HPP

#include "abs_linearize.hpp"
#include "fission_del_linearize.hpp"
#include "classical_linearize.hpp"

namespace algo { 

namespace linearize { 

template<class graph_pack_t>
struct GeneralLinearize : public algo::linearize::AbsLinearize<graph_pack_t> {
  using mcolor_t = typename graph_pack_t::mcolor_t;
  using edge_t = typename graph_pack_t::edge_t;
  using twobreak_t = typename graph_pack_t::twobreak_t; 
  using transform_t = typename graph_pack_t::transform_t;
  using partgraph_t = typename graph_pack_t::partgraph_t;   
  using genome_t = structure::Genome;
  using chromosome_t = structure::Chromosome; 
  using change_history_t = std::pair<transform_t, transform_t>;

  explicit GeneralLinearize(graph_pack_t const & gp) 
  : AbsLinearize<graph_pack_t>(gp)
  , deletion_linearize(gp)
  , classical_linearize(gp)
  {
  }

  transform_t linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q, size_t diff_chromosomes) const override; 
  
private: 
  /**
  * Split input transformation on deletion transformation -> classical transformation -> insertion transformation. 
  * Return tuple with deletions, two-breaks, insertions operations. 
  */
  transform_t get_safely_deletions(transform_t & transformation) const;
  transform_t get_safely_insertions(transform_t & transformation) const;

  void move_transformation_across_transformation(transform_t & target, transform_t & across, transform_t & source) const;

private:
	IndelLinearize<graph_pack_t> deletion_linearize;  
	ClassicalLinearize<graph_pack_t> classical_linearize; 
}; 

template<class graph_pack_t>
typename GeneralLinearize<graph_pack_t>::transform_t GeneralLinearize<graph_pack_t>::get_safely_deletions(transform_t & transformation) const { 
  transform_t deletions; 

  for (auto twobreak = transformation.begin(); twobreak != transformation.end();) { 
    vertex_t const & p = twobreak->get_vertex(0); vertex_t const & q = twobreak->get_vertex(1);
    vertex_t const & x = twobreak->get_vertex(2); vertex_t const & y = twobreak->get_vertex(3);
    
    if ((p != Infty && x != Infty && p == this->graph_pack.graph.get_obverse_vertex(x)) || 
        (q != Infty && y != Infty && q == this->graph_pack.graph.get_obverse_vertex(y))) { 
      auto second_j = twobreak; 

      if (twobreak != transformation.begin()) {  
        for (--twobreak; twobreak != transformation.begin(); --twobreak) {
          deletion_linearize.move_deletion_to_begin(transformation, twobreak, second_j);
          second_j = twobreak; 
        }
        deletion_linearize.move_deletion_to_begin(transformation, transformation.begin(), second_j);
      } 

      deletions.push_back(transformation.front());    
      transformation.pop_front(); 
      twobreak = transformation.begin(); 
    } else { 
      ++twobreak;
    }
  }

  return deletions;
}

template<class graph_pack_t>
typename GeneralLinearize<graph_pack_t>::transform_t GeneralLinearize<graph_pack_t>::get_safely_insertions(transform_t & transformation) const { 
  transform_t insertions;

  for (auto twobreak = transformation.begin(); twobreak != transformation.end();) { 
    vertex_t const & p = twobreak->get_vertex(0); vertex_t const & q = twobreak->get_vertex(1);
    vertex_t const & x = twobreak->get_vertex(2); vertex_t const & y = twobreak->get_vertex(3);

    if ((p != Infty && q != Infty && p == this->graph_pack.graph.get_obverse_vertex(q)) || 
        (x != Infty && y != Infty && x == this->graph_pack.graph.get_obverse_vertex(y))) { 
      auto first_j = twobreak; 
      if (twobreak != (--transformation.end())) {  
        for (++twobreak; twobreak != transformation.end(); ++twobreak) {
          deletion_linearize.move_insertion_to_end(transformation, first_j, twobreak);
          first_j = twobreak; 
        }
      } 
      insertions.push_front(transformation.back());
      transformation.pop_back(); 
      twobreak = transformation.begin(); 
    } else { 
      ++twobreak;
    }
  } 

  return insertions;
}

template<class graph_pack_t>
void GeneralLinearize<graph_pack_t>::move_transformation_across_transformation(transform_t & target, transform_t & across, transform_t & source) const { 
  while (!source.empty()) { 
    across.push_back(source.front());
    source.pop_front();

    auto first_j = (--across.end()); 
    if (first_j != across.begin()) {
      auto second_j = first_j;
      for (--first_j; first_j != across.begin(); --first_j) {
        this->swap_two_twobreaks(across, first_j, second_j);
        second_j = first_j; 
      } 
      this->swap_two_twobreaks(across, first_j, second_j);
    } 

    target.push_back(across.front()); 
    across.pop_front(); 
  }
}

template<class graph_pack_t>
typename graph_pack_t::transform_t GeneralLinearize<graph_pack_t>::linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q, size_t diff_chromosomes) const { 
	using namespace linearize; 

 	transform_t del_transform = get_safely_deletions(transform);
  transform_t ins_transform = get_safely_insertions(transform); 

  //Go to algorithm
  size_t c_P = count_circular_chromosome(this->graph_pack, P);
  partgraph_t PP = apply_transformation(P, del_transform); 
  size_t c_PP = count_circular_chromosome(this->graph_pack, PP);     
  partgraph_t QQ = apply_transformation(PP, transform); 
  size_t c_QQ = count_circular_chromosome(this->graph_pack, QQ); 
  assert(apply_transformation(QQ, ins_transform) == Q);
  size_t c_Q = count_circular_chromosome(this->graph_pack, Q);
  assert(c_P > c_Q); 
  INFO("Number of circular chromosomes" << c_P << " " << c_PP << " " << c_QQ << " " << c_Q)

  //Run classical linearization.
  transform_t linearize_transform; 
  if (c_PP > c_QQ) {
    INFO("Start linearize algortihm on classic history")
    linearize_transform = classical_linearize.linearize(PP, transform, QQ, diff_chromosomes);
    assert(linearize_transform.size() == std::min(c_PP - c_QQ, diff_chromosomes)); 
    diff_chromosomes -= std::min(c_PP - c_QQ, diff_chromosomes);
    INFO("Finish linearize algortihm on classic history")
    
  } 

  //Run deletion linearization
  transform_t del_linearize_transform; 
  if (c_P > c_PP) {
    INFO("Start linearize algortihm on deletion history")
    del_linearize_transform = deletion_linearize.linearize(P, del_transform, PP, diff_chromosomes);
    assert(del_linearize_transform.size() == std::min(c_P - c_PP, diff_chromosomes)); 
    diff_chromosomes -= std::min(c_P - c_PP, diff_chromosomes);
    INFO("Finish linearize algortihm on deletion history")
  }

  move_transformation_across_transformation(del_linearize_transform, del_transform, linearize_transform);
  assert(del_linearize_transform.size() == (c_P - c_Q)); 
  transform.splice(transform.end(), ins_transform); transform.splice(transform.begin(), del_transform);
  
  return del_linearize_transform;
} 

} 

} 


#endif