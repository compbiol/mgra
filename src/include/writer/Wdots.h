#ifndef WDOTS_H_
#define WDOTS_H_

/*
 * Fixme remove RGDCoeff graph picture
 *
 */
namespace writer {
  template<class graph_t, class conf_t>
  struct Wdots { 
    typedef typename graph_t::edge_t edge_t;
    typedef typename graph_t::mcolor_type mcolor_t;
    typedef typename graph_t::mularcs_t mularcs_t;

    Wdots() 
    : m_debug(false)
    {
    }

    void init(std::string const & path, std::string const & graphname, bool debug) {
      m_path = path; 
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
        save_dot(graph, path::append_path(m_path, "debug"), dotname);    
      }
    }

    void save_components(graph_t const & graph, size_t stage);
    void write_legend_dot();

    void save_median(std::string const & name, std::vector<std::pair<edge_t, mcolor_t> > const & edges);

  private: 
    void save_dot(graph_t const & graph, std::string const & path, std::string const & dotname);

  private: 
    std::string m_path;
    std::string m_colorscheme;
    std::string m_graphname;
    bool m_debug;
  }; 
} 

template<class graph_t, class conf_t>
void writer::Wdots<graph_t, conf_t>::save_median(std::string const & dotname, std::vector<std::pair<edge_t, mcolor_t> > const & edges) { 
  std::ofstream dot(path::append_path(m_path, dotname));
  
  dot << "graph {" << std::endl;
  
  if (!cfg::get().colorscheme.empty()) { 
    dot << "edge [colorscheme=" << cfg::get().colorscheme << "];" << std::endl;
  } 

  int infv = 0;
  for(auto const & edge : edges) { 
    mcolor_t const & color = edge.second;
    for(auto ic = color.cbegin(); ic != color.cend(); ++ic) {
      for (size_t i = 0; i < ic->second; ++i) { 
      
        dot << "\t\"" << edge.first.first << "\"\t--\t\"";
        if (edge.first.second == Infty) {
          if (ic == color.cbegin()) { 
            --infv;
          } 
          dot << infv << "\"\t[len=0.75,";
        } else { 
          dot << edge.first.second << "\"\t[";
        } 
        dot << "color=" <<  cfg::get().get_RGBcolor(ic->first) << ", penwidth=3];" << std::endl;
      }
    } 
  }  

  for(int i = infv; i < 0; ++i) {
    dot << "\t\"" << i << "\"\t[shape=point,color=black];" << std::endl;
  }
  
  dot << "}" << std::endl;
  dot.close();
}

template<class graph_t, class conf_t>
void writer::Wdots<graph_t, conf_t>::save_dot(graph_t const & graph, std::string const & path, std::string const & dotname) { 
  std::ofstream dot(path::append_path(path, dotname));

  dot << "graph {" << std::endl;
  if (!cfg::get().colorscheme.empty()) { 
    dot << "edge [colorscheme=" << cfg::get().colorscheme << "];" << std::endl;
  } 

  int infv = 0;
  std::unordered_set<vertex_t> mark; // vertex set
  for(auto const & x : graph) { 
    mularcs_t const & Mx = graph.get_all_adjacent_multiedges(x);

    if (Mx.size() == 1 && Mx.union_multicolors() == graph.get_complete_color()) { 
      continue; // trivial cycle
    } 

    for(auto im = Mx.cbegin(); im != Mx.cend(); ++im) {
      vertex_t const & y = im->first;
      if (mark.find(y) != mark.end()) { 
        continue; // already output
      }    

      //std::cerr << x << " " << is_pseudo << " " << y << std::endl;

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

            dot << "color=" <<  cfg::get().get_RGBcolor(ic->first) << ", penwidth=3];" << std::endl;
          }
        }
      } else { 
        auto split_colors = graph.split_color(im->second);

        if (split_colors.size() == 2) { 
          auto const & C = *split_colors.begin();
          for(auto ic = C.cbegin(); ic != C.cend(); ++ic) {
            for (size_t i = 0; i < ic->second; ++i) { 
              //************** output edge (x,y) **************** 
              dot << "\t\"" << x << "\"\t--\t\"";
              if (y == Infty) {
                if (ic == C.cbegin()) { 
                  --infv;
                } 
                dot << infv << "\"\t[len=0.75,";
              } else { 
                dot << y << "\"\t[";
              }

              dot << "color=" <<  cfg::get().get_RGBcolor(ic->first) << ", style=dashed];" << std::endl;   

            }
          } 

          
          auto const & C1 = *(++split_colors.cbegin());
          for(auto ic = C1.cbegin(); ic != C1.cend(); ++ic) {
            for (size_t i = 0; i < ic->second; ++i) { 
              //************** output edge (x,y) **************** 
              dot << "\t\"" << x << "\"\t--\t\"";
              if (y == Infty) {
                if (ic == C1.cbegin()) { 
                  --infv;
                } 
                dot << infv << "\"\t[len=0.75,";
              } else { 
                dot << y << "\"\t[";
              }

              dot << "color=" <<  cfg::get().get_RGBcolor(ic->first) << "];" << std::endl;   
            }
          }           
        } else { 
          mcolor_t const & C = im->second;
          size_t number_splits = graph.split_color(im->second).size();
          bool out = false;
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

              if (!out && number_splits != 1) {
                out = true; 

                dot << "color=" <<  cfg::get().get_RGBcolor(ic->first) << ", label=" << number_splits << "];" << std::endl;  
              } else { 
                dot << "color=" <<  cfg::get().get_RGBcolor(ic->first) << "];" << std::endl;  
              } 
            } 
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
    std::ofstream dot(path::append_path(m_path, namefile));

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

      mularcs_t const & Mx = graph.get_all_adjacent_multiedges(x);

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
	      dot << "color=" <<  cfg::get().get_RGBcolor(ic->first) << ", penwidth=3];" << std::endl;
	    } else {
	      dot << "color=" <<  cfg::get().get_RGBcolor(ic->first) << "];" << std::endl;	
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
  std::ofstream output(path::append_path(m_path, "legend.dot"));

  output << "digraph legend {" << std::endl;
  output << "\tnode [style=filled"; 
  if (!m_colorscheme.empty()) {
    output << ", colorscheme=" << m_colorscheme;
  } 
  output << "];" << std::endl;

  for (size_t j = 0; j < cfg::get().get_count_genomes(); ++j) {
    output << "\t\"" << cfg::get().get_priority_name(j) << "\" [fillcolor=" 
      << cfg::get().get_RGBcolor(j)  << "];" << std::endl;
  } 

  std::vector<std::string> info;
  for(auto const & tree : cfg::get().phylotrees) {
    tree.get_nodes(info);
  }
 
  for(auto const & str : info) {
    output << str << std::endl;
  } 

  output << "}" << std::endl;
  output.close();
} 

#endif

