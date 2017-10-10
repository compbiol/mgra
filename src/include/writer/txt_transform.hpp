#ifndef TXT_TRANSFORMATION_HPP
#define TXT_TRANSFORMATION_HPP

namespace writer {

template <class graph_pack_t>
struct TXT_transformation {
	using mcolor_t = typename graph_pack_t::mcolor_t;
	using transform_t = typename graph_pack_t::transform_t;
	using partgraph_t = typename graph_pack_t::partgraph_t;

  TXT_transformation(std::string const & path, graph_pack_t const & gp) 
  : m_path(path)
  , graph_pack(gp)
  { 
  }

  void save_full_colors_transformation(transform_t const & transform) { 
    std::string namefile = "full_history.txt";
    std::string new_path = path::append_path(m_path, namefile);
    std::ofstream out(new_path);
    ;//write_transformation(out, transform); 
    out.close();  
  } 
  
  void save_transformation(std::pair<mcolor_t, mcolor_t> const & branch, transform_t const & transform) { 
  	std::string namefile = cfg::get().mcolor_to_name(branch.first) + "--" + cfg::get().mcolor_to_name(branch.second) + ".trs";
  	std::string new_path = path::append_path(m_path, path::append_path("transformations", namefile));
  	std::ofstream out(new_path);
  	write_transformation(out, transform); 
  	out.close();  
  }

  void save_reverse_transformation(std::pair<mcolor_t, mcolor_t> const & branch, transform_t const & transform) { 
  	std::string namefile = cfg::get().mcolor_to_name(branch.second) + "--" + cfg::get().mcolor_to_name(branch.first) + ".trs";
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
    	  if (p == graph_pack.graph.get_obverse_vertex(x) && graph_pack.is_prosthetic_chromosome(p, x)) { 
    	    out << " # deletion"; 
        } else if (q == graph_pack.graph.get_obverse_vertex(y) && graph_pack.is_prosthetic_chromosome(q, y)) {
    	    out << " # deletion"; 
    	  } else if (p == graph_pack.graph.get_obverse_vertex(q) && graph_pack.is_prosthetic_chromosome(p, q)) { 
    	    out << " # insertion"; 
        } else if (x == graph_pack.graph.get_obverse_vertex(y) && graph_pack.is_prosthetic_chromosome(x, y)) {
    	    out << " # insertion"; 
    	  } 
    	}
    	out << std::endl;  
    }
	}
private: 
	std::string m_path;
	graph_pack_t const & graph_pack;
}; 	

} 

#endif
