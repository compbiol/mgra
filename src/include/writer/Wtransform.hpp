#ifndef WTRANSFORMATION_HPP
#define WTRANSFORMATION_HPP

namespace writer {

template <class graph_t>
struct Wtransformation {
	using mcolor_t = typename graph_t::mcolor_type;
	using transform_t = typename graph_t::transform_t;
	using partgraph_t = typename graph_t::partgraph_t;

  Wtransformation(std::string const & path, graph_t const & graph) 
  : m_path(path)
  , m_graph(graph)
  , bad_edges(graph.get_bad_edges())
  { 
  }

  void save_transformation(std::pair<mcolor_t, mcolor_t> const & branch, transform_t const & transform) { 
  	std::string namefile = cfg::get().mcolor_to_name(branch.first) + "--" 
  													+ cfg::get().mcolor_to_name(branch.second) + ".trs";
  	std::string new_path = path::append_path(m_path, path::append_path("transformations", namefile));
  	std::ofstream out(new_path);
  	write_transformation(out, transform); 
  	out.close();  
  }

  void save_reverse_transformation(std::pair<mcolor_t, mcolor_t> const & branch, transform_t const & transform) { 
  	std::string namefile = cfg::get().mcolor_to_name(branch.second) + "--" 
  													+ cfg::get().mcolor_to_name(branch.first) + ".trs";
  	std::string new_path = path::append_path(m_path, path::append_path("transformations", namefile));
  	transform_t reverse_transform;
   	for (auto twobreak = transform.rbegin(); twobreak != transform.rend(); ++twobreak) {
   		reverse_transform.push_back(twobreak->inverse()); 
   	}
   	std::ofstream out(new_path);
  	write_transformation(out, reverse_transform); 
  	out.close();  
  }

private: 
	void write_transformation(std::ofstream & out, transform_t const & transform) { 
		for(auto const & event : transform) {
      vertex_t const & p = event.get_vertex(0);
      vertex_t const & q = event.get_vertex(1);
      vertex_t const & x = event.get_vertex(2);
      vertex_t const & y = event.get_vertex(3);
    
    	out << "(" << p << ", " << q << ") x (" << x << ", " << y << ") " << cfg::get().mcolor_to_name(event.get_mcolor()); 

    	if (p != Infty && q != Infty && x != Infty && y != Infty) { 
    	  if (p == m_graph.graph.get_obverse_vertex(x) && bad_edges.defined(p, x)) { 
    	    out << " # deletion"; 
        } else if (q == m_graph.graph.get_obverse_vertex(y) && bad_edges.defined(q, y)) {
    	    out << " # deletion"; 
    	  } else if (p == m_graph.graph.get_obverse_vertex(q) && bad_edges.defined(p, q)) { 
    	    out << " # insertion"; 
        } else if (x == m_graph.graph.get_obverse_vertex(y) && bad_edges.defined(x, y)) {
    	    out << " # insertion"; 
    	  } 
    	}
    	out << std::endl;  
    }
	}
private: 
	std::string m_path;
	graph_t const & m_graph;
	partgraph_t bad_edges;
}; 	

} 

#endif