#ifndef WDOTS_H_
#define WDOTS_H_

namespace writer {
  template<class graph_t, class conf_t>
  struct Wdots { 
    typedef structure::Mcolor mcolor_t;
    typedef structure::Mularcs<mcolor_t> mularcs_t;

    Wdots(conf_t const & cfg) 
    : m_cfg(cfg)
    , m_debug(false)
    {
    }

    void init(fs::path const & path, std::string const & colorscheme, std::string const & graphname, bool debug) {
      m_path = path; 
      m_colorscheme = colorscheme; 
      m_graphname = graphname;
      m_debug = debug;
    }     

    // Save .dot file and output statistics of synteny blocks representing breakpoints
    void save_final_dot(graph_t const & graph) {
      save_dot(graph, m_path, "last_graph.dot");
    }

    void save_dot(graph_t const & graph, size_t stage) {
      if (m_debug) {
        std::string dotname = m_graphname + std::to_string(stage) + ".dot";
        save_dot(graph, m_path / "debug", dotname);    
      }
    }

    void save_components(graph_t const & graph, size_t stage);
    void write_legend_dot();

  private: 
    void save_dot(graph_t const & graph, fs::path const & path, std::string const & dotname);
  private: 
    conf_t const & m_cfg;
    fs::path m_path;
    std::string m_colorscheme;
    std::string m_graphname;
    bool m_debug;
  }; 
} 

template<class graph_t, class conf_t>
void writer::Wdots<graph_t, conf_t>::save_dot(graph_t const & graph, fs::path const & path, std::string const & dotname) { 
  fs::ofstream dot(path / dotname);

  dot << "graph {" << std::endl;
  if (!m_colorscheme.empty()) { 
    dot << "edge [colorscheme=" << m_colorscheme << "];" << std::endl;
  } 

  int infv = 0;
  std::unordered_set<vertex_t> mark; // vertex set
  for(auto const & x : graph) { 
    mularcs_t const & Mx = graph.get_adjacent_multiedges(x);

    if (Mx.size() == 1 && Mx.union_multicolors() == graph.get_complete_color()) { 
      continue; // trivial cycle
    } 

    for(auto im = Mx.cbegin(); im != Mx.cend(); ++im) {
      vertex_t const & y = im->first;

      if (mark.find(y) != mark.end()) { 
        continue; // already output
      }    

      if (graph.is_vec_T_consistent_color(im->second)) { 
        mcolor_t const & C = im->second;
        for(auto ic = C.cbegin(); ic != C.cend(); ++ic) {
          for (size_t i = 0; i < ic->second; ++i) { 
            /*************** output edge (x,y) **************** */
            dot << "\t\"" << x << "\"\t--\t\"";
            if (y == Infty) {
              if (ic == C.cbegin()) { 
                --infv;
              } 
              dot << infv << "\"\t[len=0.75,";
            } else { 
              dot << y << "\"\t[";
            } 

            dot << "color=" <<  m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * (ic->first)) << ", penwidth=3];" << std::endl;
          }
        }
      } else { 
        auto split_colors = graph.split_color(im->second);
        if (split_colors.size() == 2) { 
          auto const & C = *split_colors.begin();
          for(auto ic = C.cbegin(); ic != C.cend(); ++ic) {
            for (size_t i = 0; i < ic->second; ++i) { 
              /*************** output edge (x,y) **************** */
              dot << "\t\"" << x << "\"\t--\t\"";
              if (y == Infty) {
                if (ic == C.cbegin()) { 
                  --infv;
                } 
                dot << infv << "\"\t[len=0.75,";
              } else { 
                dot << y << "\"\t[";
              }

              dot << "color=" <<  m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * (ic->first)) << ", style=dashed];" << std::endl;   
            }
          } 

          
          auto const & C1 = *(++split_colors.cbegin());
          for(auto ic = C1.cbegin(); ic != C1.cend(); ++ic) {
            for (size_t i = 0; i < ic->second; ++i) { 
              /*************** output edge (x,y) **************** */
              dot << "\t\"" << x << "\"\t--\t\"";
              if (y == Infty) {
                if (ic == C1.cbegin()) { 
                  --infv;
                } 
                dot << infv << "\"\t[len=0.75,";
              } else { 
                dot << y << "\"\t[";
              }

              dot << "color=" <<  m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * (ic->first)) << "];" << std::endl;   
            }
          }           
        } else { 
          mcolor_t const & C = im->second;
          for(auto ic = C.cbegin(); ic != C.cend(); ++ic) {
            for (size_t i = 0; i < ic->second; ++i) { 
              /*************** output edge (x,y) **************** */
              dot << "\t\"" << x << "\"\t--\t\"";
              if (y == Infty) {
                if (ic == C.cbegin()) { 
                  --infv;
                } 
                dot << infv << "\"\t[len=0.75,";
              } else { 
                dot << y << "\"\t[";
              } 

              dot << "color=" <<  m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * (ic->first)) << "];" << std::endl;  
            } 
          }
        }
      } 
      /*for(auto ic = C.cbegin(); ic != C.cend(); ++ic) {
        for (size_t i = 0; i < ic->second; ++i) { 
          /*************** output edge (x,y) **************** 
          dot << "\t\"" << x << "\"\t--\t\"";
          if (y == Infty) {
            if (ic == C.cbegin()) { 
              --infv;
            } 
            dot << infv << "\"\t[len=0.75,";
          } else { 
            dot << y << "\"\t[";
          } 

          if (vec_T_color) {
            dot << "color=" <<  m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * (ic->first)) << ", penwidth=3];" << std::endl;
          } else {
            dot << "color=" <<  m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * (ic->first)) << "];" << std::endl;	
          }
        } 
      }*/
    }
    mark.insert(x);
  }

  for(int i = infv; i < 0; ++i) {
    dot << "\t\"" << i << "\"\t[shape=point,color=black];" << std::endl;
  }

  dot << "}" << std::endl;
  dot.close();
} 

template<class graph_t, class conf_t>
void writer::Wdots<graph_t, conf_t>::save_components(graph_t const & graph, size_t stage) { 
  std::string dotname = m_graphname + std::to_string(stage);
  std::map<vertex_t, std::set<vertex_t> > components = graph.split_on_components().get_eclasses(); 
  
  size_t i = 0; 
  for(auto it = components.cbegin(); it != components.cend(); ++it) { 
    auto const & current = it->second;
    if (current.size() <= 2) {
      continue;
    }
  
    std::string namefile = dotname + "_" + std::to_string(++i) + ".dot"; 
    fs::ofstream dot(m_path / namefile);

    dot << "graph {" << std::endl;
    if (!m_colorscheme.empty()) { 
      dot << "edge [colorscheme=" << m_colorscheme << "];" << std::endl;
    } 

    int infv = 0;
    std::unordered_set<vertex_t> mark; // vertex set
    for(auto is = current.cbegin(); is != current.cend(); ++is) {
      vertex_t const & x = *is;
  
      if (x == Infty) { 
      	continue;
      } 

      mularcs_t const & Mx = graph.get_adjacent_multiedges(x);

      if (Mx.union_multicolors() == graph.get_complete_color()) { 
        continue; // trivial cycle
      } 
      
      for(auto im = Mx.cbegin(); im != Mx.cend(); ++im) {
	vertex_t const & y = im->first;
	
	if (mark.find(y) != mark.end()) { 
	  continue; // already output
	} 

	mcolor_t const & C = im->second;
	bool vec_T_color = graph.is_vec_T_consistent_color(C);
	for(auto ic = C.cbegin(); ic != C.cend(); ++ic) {
	  for (size_t i = 0; i < ic->second; ++i) { 
	    dot << "\t\"" << x << "\"\t--\t\"";
	    if (y == Infty) {
	      if (ic == C.cbegin()) { 
		      --infv;
	      } 
	      dot << infv << "\"\t[len=0.75,";
	    } else { 
	      dot << y << "\"\t[";
	    } 
	    if (vec_T_color) {
	      dot << "color=" <<  m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * (ic->first)) << ", penwidth=3];" << std::endl;
	    } else {
	      dot << "color=" <<  m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * (ic->first)) << "];" << std::endl;	
	    }
	  }
	}
      }
      mark.insert(x);
    }

    for(int i = infv; i < 0; ++i) {
      dot << "\t\"" << i << "\"\t[shape=point,color=black];" << std::endl;
    }

    dot << "}" << std::endl;
    dot.close();
  } 
} 

template<class graph_t, class conf_t>
void writer::Wdots<graph_t, conf_t>::write_legend_dot() { 
  fs::ofstream output(m_path / "legend.dot");

  output << "digraph legend {" << std::endl;
  output << "\tnode [style=filled"; 
  if (!m_colorscheme.empty()) {
    output << ", colorscheme=" << m_colorscheme;
  } 
  output << "];" << std::endl;

  for (size_t j = 0; j < m_cfg.get_count_genomes(); ++j) {
    output << "\t\"" << m_cfg.get_priority_name(j) << "\" [fillcolor=" 
      << m_cfg.get_RGBcolor(m_cfg.get_RGBcoeff() * j)  << "];" << std::endl;
  } 

  std::vector<std::string> info;
  for(auto it = m_cfg.cbegin_trees(); it != m_cfg.cend_trees(); ++it) {
    it->get_nodes(info);
  }
 
  for(auto const & str : info) {
    output << str << std::endl;
  } 

  output << "}" << std::endl;
  output.close();
} 

#endif

