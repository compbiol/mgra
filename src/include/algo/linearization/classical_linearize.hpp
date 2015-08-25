#ifndef CLASSICAL_LINEARIZE_HPP
#define CLASSICAL_LINEARIZE_HPP

namespace algo { 

namespace linearize { 

template<class graph_pack_t>
struct ClassicalLinearize : public algo::linearize::AbsLinearize<graph_pack_t> {
	using mcolor_t = typename graph_pack_t::mcolor_type;
  using edge_t = typename graph_pack_t::edge_t;
  using twobreak_t = typename graph_pack_t::twobreak_t; 
  using transform_t = typename graph_pack_t::transform_t;
  using partgraph_t = typename graph_pack_t::partgraph_t; 
  using citer_transform = typename transform_t::iterator; 
  using dependence_type = typename twobreak_t::dependence_type;
  using genome_t = structure::Genome;
  using chromosome_t = structure::Chromosome; 
  using change_history_t = std::pair<transform_t, transform_t>;

  explicit ClassicalLinearize(graph_pack_t const & gp) 
  : AbsLinearize<graph_pack_t>(gp)
  {
  }

  change_history_t linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q) const override; 

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

	/**
	 * Start subranges [....][<....][<.....][<.....][.....>]
	 */
	std::list<citer_transform> split_range(citer_transform start, partgraph_t current, citer_transform finish) const {
		size_t c_q0 = count_circular_chromosome(this->graph_pack, current);	size_t c_q1 = c_q0;

		std::list<citer_transform> start_subranges; 
		for (auto iter = start; iter != finish; ++iter) { 
			iter->apply_single(current);
			c_q0 = c_q1; c_q1 = count_circular_chromosome(this->graph_pack, current); 
			if (c_q0 < c_q1) { 
				start_subranges.push_back(iter);
			}
		}

		return start_subranges;
	} 

	citer_transform move_down(citer_transform range_start, partgraph_t current, citer_transform range_finish, citer_transform finish) const { 
		size_t c_q0 = count_circular_chromosome(this->graph_pack, current);	size_t c_q1 = c_q0;
		partgraph_t P = current;

		for (auto iter = range_start; iter != range_finish && iter != finish;) { 
			iter->apply_single(P);
			c_q0 = c_q1; c_q1 = count_circular_chromosome(this->graph_pack, P);

			if (c_q0 > c_q1) { 
				auto first = iter; auto second = iter; 

				while (first != range_start) {
					second = first; --first; 
					finish = this->swap_two_twobreaks(first, second, finish);
				}

				range_start->apply_single(current);
				iter = (++range_start);
				P = current;
			} else { 
				++iter;
			}
		}

		return finish;
	}

	std::tuple<citer_transform, citer_transform> move_up(citer_transform start_range, partgraph_t Q, citer_transform finish_range, citer_transform finish) const { 
		auto first = start_range;	auto second = ++start_range; auto third = ++start_range;
		
		partgraph_t P = Q;
		first->apply_single(Q);	second->apply_single(Q);
		size_t c_q0 = count_circular_chromosome(this->graph_pack, Q);
		third->apply_single(Q);
		size_t c_q1 = count_circular_chromosome(this->graph_pack, Q);
		
		assert(c_q0 > c_q1);
		while ((third != finish_range) && (c_q0 > c_q1)) {
			finish = this->swap_three_twobreaks(P, first, second, third, finish);
			
			first->apply_single(P);
			first = second; second = third; ++third; 

			if (third != finish) { 
				third->apply_single(Q);
				c_q0 = c_q1; c_q1 = count_circular_chromosome(this->graph_pack, Q);
			} 
		}

		return std::make_pair(first, finish);
	} 

	citer_transform process_range_by_range(citer_transform start_range, partgraph_t current, citer_transform finish_range, citer_transform finish,
										std::list<citer_transform> const & small_ranges) const { 
		auto finish_small_range = finish_range;
		for (auto local_start_range = small_ranges.crbegin(); local_start_range != small_ranges.crend(); ++local_start_range) { 
			partgraph_t P = current; 
			for (auto iter = start_range; (iter != finish_range) && (iter != *local_start_range); ++iter) { 
				iter->apply_single(P);
			} 

			std::cerr << "Start new range " << std::endl;			
			partgraph_t T = P;
			auto temp = *local_start_range;
			std::cerr << count_circular_chromosome(this->graph_pack, T); 
			for (auto it = temp; it != finish_small_range; ++it) { 
				it->apply_single(T);
				std::cerr << "->" << count_circular_chromosome(this->graph_pack, T); 
			}
			std::cerr << std::endl;
			
			// move to begin all >
			auto start_small_range = *local_start_range;
			(*local_start_range)->apply_single(P); 
			finish = move_down((++start_small_range), P, finish_small_range, finish);

			std::cerr << "Result after move down" << std::endl;
			T = P; 
			(*local_start_range)->inverse().apply_single(T); 
			std::cerr << count_circular_chromosome(this->graph_pack, T); 
			for (auto it = *local_start_range; it != finish_small_range; ++it) { 
				it->apply_single(T);
				std::cerr << "->" << count_circular_chromosome(this->graph_pack, T);
			}
			std::cerr << std::endl;
			 
			// swap all <>>
			std::cerr << "Start move up all <>> " << std::endl;	
			start_small_range = *local_start_range; 
			(*local_start_range)->inverse().apply_single(P); 
			std::tie(finish_small_range, finish) = move_up(start_small_range, P, finish_small_range, finish);
			std::cerr << "Finish move up all <>> " << std::endl;	

			std::cerr << "Result after move up" << std::endl;
			T = P; 
			std::cerr << count_circular_chromosome(this->graph_pack, T); 
			for (auto it = *local_start_range; it != finish_small_range; ++it) { 
				it->apply_single(T);
				std::cerr << "->" << count_circular_chromosome(this->graph_pack, T);
			}
			std::cerr << std::endl;
		}		

		if (*small_ranges.begin() != start_range) { 
			finish = move_down(start_range, current, finish_small_range, finish);
		}

		return finish;
	}

	citer_transform one_step_induction(citer_transform start, partgraph_t current, citer_transform finish) const; 	

	citer_transform swap_three_twobreaks(partgraph_t const & current, citer_transform first, citer_transform second, citer_transform third, citer_transform finish) const { 		
		std::cerr << "First break: " << first->get_vertex(0) << " " << first->get_vertex(1) << " " 
				<< first->get_vertex(2) << " " << first->get_vertex(3) << std::endl;

		std::cerr << "Second break: " << second->get_vertex(0) << " " << second->get_vertex(1) << " " 
				<< second->get_vertex(2) << " " << second->get_vertex(3) << std::endl;

		std::cerr << "Third break: " << third->get_vertex(0) << " " << third->get_vertex(1) << " " 
				<< third->get_vertex(2) << " " << third->get_vertex(3) << std::endl;
	  
	  if (first->is_dependent(*second) == twobreak_t::independent) { 
	  	std::cerr << "SWAP THREE: first and second independent case "; 
	  	if (is_belong_one_chromosome_indep(current, *second)) { 
	  		std::cerr << "and first and second independent case and belong one chromosome" << std::endl; 
	    	//Theorem 4. if a and b belong to one chromosome.
	    	//Theorem 4. if a and b belong to one linear chromosome.
	      finish = this->swap_two_twobreaks(second, third, finish);
	      finish = this->swap_two_twobreaks(first, second, finish);
	    } else {
	    	std::cerr << "and first and second independent case: basic swap" << std::endl; 
	     	//Theorem 4. if a and b belong to two different chromosomes.
	     	std::iter_swap(first, second);
	    }
	  } else if (first->is_dependent(*second) == twobreak_t::weakly_dependent) {
	  	//std::cerr << "SWAP THREE: first and second weakly dependent case "; 	  	
	    if (is_belong_one_chromosome_dep(current, *first, *second)) { 
	    	//std::cerr << "and first and second weakly dependent case: linear or belong one" << std::endl; 
	    	//Theorem 5. if a and b and d belong to one chromosome. 
	    	//Theorem 5. if a and b and d belong to two different chromosomes and both linear.  
	      finish = this->swap_two_twobreaks(second, third, finish);
	      finish = this->swap_two_twobreaks(first, second, finish);
	    } else {  
	    	//std::cerr << "and first and second weakly dependent case: basic swap" << std::endl; 
	      //Theorem 5. if a and b and d belong to two different chromosome and one of circular.
	      finish = this->swap_two_twobreaks(first, second, finish);
	    } 
	  } else { 
	  	finish = this->swap_two_twobreaks(first, second, finish);
	  } 
	  return finish;
	}

	bool is_belong_one_chromosome_indep(partgraph_t P, twobreak_t const & twobreak) const; // const &
	bool is_belong_one_chromosome_dep(partgraph_t const & P, twobreak_t const & first, twobreak_t const & second) const;
}; 

template<class graph_pack_t>
bool ClassicalLinearize<graph_pack_t>::is_belong_one_chromosome_indep(partgraph_t P, twobreak_t const & twobreak) const { //const & 
	auto const & a = twobreak.get_arc(0); auto const & b = twobreak.get_arc(1);
	std::unordered_set<vertex_t> chromosome_set; 
	auto chr1 = get_chromosome_by_edge(this->graph_pack, P, a, chromosome_set);

	/*std::cerr << "Test in indep belong " << std::endl;
	std::unordered_set<vertex_t> chromosome_set1; 
	auto chr3 = get_chromosome_by_edge(this->graph_pack, P, b, chromosome_set1);
	std::cerr << chr1.is_circular() << " " << chr3.is_circular() << std::endl;
	std::cerr << (chromosome_set1 == chromosome_set) << std::endl;
	*/

	assert(a.second == Infty || chromosome_set.count(a.second) != 0);
	if ((b.first == Infty || chromosome_set.count(b.first) != 0) 
			&& (b.second == Infty || chromosome_set.count(b.second) != 0) 
			&& P.defined(b.first, b.second)) { 
		//std::cerr << "we are here " << std::endl;
		return true; 
	}

	auto chr2 = get_chromosome_by_edge(this->graph_pack, P, b, chromosome_set);
	return !(chr1.is_circular() || chr2.is_circular());	
}

template<class graph_pack_t>
bool ClassicalLinearize<graph_pack_t>::is_belong_one_chromosome_dep(partgraph_t const & P, twobreak_t const & first, twobreak_t const & second) const { 
	auto check_lambda = [&] (size_t ind1, size_t ind2, size_t ind3) -> bool {
  	if (second.get_arc(ind1) != edge_t(Infty, Infty)) { 
    	if (second.get_arc(ind1) == std::make_pair(first.get_vertex(ind2), first.get_vertex(ind3))) { 
      	return true;
    	} 
    	if (second.get_arc(ind1) == std::make_pair(first.get_vertex(ind3), first.get_vertex(ind2))) { 
        return true;
	    }
  	} 
  	return false;
	};

	auto const & a = first.get_arc(0); auto const & b = first.get_arc(1);
	edge_t c; edge_t d; 

	if (check_lambda(0, 0, 2) || check_lambda(0, 1, 3)) { 
		c = second.get_arc(0); d = second.get_arc(1);
	} else if (check_lambda(1, 0, 2) || check_lambda(1, 1, 3)) { 
		c = second.get_arc(1); d = second.get_arc(0);
	}   
	
	std::unordered_set<vertex_t> chromosome_set; 
	auto chr1 = get_chromosome_by_edge(this->graph_pack, P, a, chromosome_set);
	assert((a.second == Infty || chromosome_set.count(a.second) != 0)); 
	assert((b.first == Infty || chromosome_set.count(b.first) != 0) && (b.second == Infty || chromosome_set.count(b.second) != 0));

	if ((d.first == Infty || chromosome_set.count(d.first) != 0) 
			&& (d.second == Infty || chromosome_set.count(d.second) != 0) 
			&& P.defined(d.first, d.second)) { 
		return true; 
	}

	auto chr2 = get_chromosome_by_edge(this->graph_pack, P, d, chromosome_set);
	return !(chr1.is_circular() || chr2.is_circular());
}

template<class graph_pack_t>
typename ClassicalLinearize<graph_pack_t>::citer_transform ClassicalLinearize<graph_pack_t>::one_step_induction(citer_transform start, partgraph_t current, citer_transform finish) const {  
	auto range = find_range(start, current, finish);

	partgraph_t T = current; 
	std::cerr << "Find big range" << std::endl;
	std::cerr << count_circular_chromosome(this->graph_pack, T); 
	for (auto it = range.first; it != range.second; ++it) { 
		it->apply_single(T);
		std::cerr << "->" << count_circular_chromosome(this->graph_pack, T);
	}
	std::cerr << std::endl;

	std::cerr << "Split range" << std::endl;
	auto small_ranges = split_range(range.first, current, range.second);
	std::cerr << "Size of ranges " << small_ranges.size() << std::endl;

	if (small_ranges.empty()) { 
		std::cerr << "Process begining range without increasing and without split small ranges" << std::endl;
		finish = move_down(range.first, current, range.second, finish);
	} else { 
		std::cerr << "Process small range by small range" << std::endl;
		finish = process_range_by_range(range.first, current, range.second, finish, small_ranges);
	}

	return finish; 
} 

template<class graph_pack_t>
typename ClassicalLinearize<graph_pack_t>::change_history_t ClassicalLinearize<graph_pack_t>::linearize(partgraph_t P, transform_t & transform, partgraph_t const & Q) const { 
	auto start = transform.begin(); auto finish = transform.end();

	size_t diff_chromosomes = count_circular_chromosome(this->graph_pack, P) - count_circular_chromosome(this->graph_pack, Q);
  for (size_t i = 0; i < diff_chromosomes; ++i) { 
  	std::cerr << "LALALA Start step induction " << count_circular_chromosome(this->graph_pack, P) << std::endl;     
    finish = one_step_induction(start, P, finish);
    start->apply_single(P);
    ++start;
    std::cerr << "LALALA Finish step induction " << count_circular_chromosome(this->graph_pack, P) << std::endl; 
  }

  transform.erase(finish, transform.end());

  transform_t replace_transformation; transform_t transformation = transform; 
	for (auto it = transform.begin(); it != start; ++it) { 
		replace_transformation.push_back(*transformation.begin());
    transformation.pop_front();  
  } 

  return std::make_pair(replace_transformation, transformation); 
}

}

} 

/*
std::cerr << "State after process range " << std::endl;
			T = current;
			std::cerr << count_circular_chromosome(this->graph_pack, T); 
			for (auto it = start_range; it != finish_range; ++it) { 
				it->apply_single(T);
				std::cerr << "->" << count_circular_chromosome(this->graph_pack, T);
			}
			std::cerr << std::endl;

		auto s = start_range;
		auto mem_T = Q;
std::cerr << "Interval after swap 3" << std::endl;
			auto T = mem_T;
			std::cerr << count_circular_chromosome(this->graph_pack, T); 
			for (auto it = s; it != finish_range; ++it) { 
				it->apply_single(T);
				std::cerr << "->" << count_circular_chromosome(this->graph_pack, T); 
			}
			std::cerr << std::endl;
			
	*/
#endif
