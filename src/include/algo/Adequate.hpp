#ifndef ADEQUATE_HPP
#define ADEQUATE_HPP

template<class graph_t>
struct Algorithm<graph_t>::Adequate : public Algorithm<graph_t>::Stage { 
	typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::edge_t edge_t; 

	explicit ProcessSimplePath(std::shared_ptr<graph_t> const & graph)
	: Stage(graph) 
	{
	}

	bool do_action() override;
  
  std::string get_name() override { 
    return "Process good and simple paths.";
  }

};

template<class graph_t>
bool Algorithm<graph_t>::Adequate::do_action() override { 
	for (vertex_t const & v : *this->graph) { 
	}
}
  
#endif 