#ifndef EXPEREM_STAGE_H_
#define EXPEREM_STAGE_H_

//experemental stage
template<class graph_t>
bool Algorithm<graph_t>::stage5_4() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 

  do {
      number_rear = 0;
        utility::equivalence<vertex_t> classes = graph->split_on_components();
        std::map<vertex_t, std::set<vertex_t> > components = classes.get_eclasses<std::set<vertex_t> >();


      for (auto const & component : components) {
        std::set<vertex_t> const & verteces = component.second;
        std::unordered_set<vertex_t> processed;

        std::map<mcolor_t, std::vector<arc_t> > how_many_edges;
        for (auto const & vertex : verteces) { 
          mularcs_t mularcs = graph->get_adjacent_multiedges_with_info(vertex);
          for (auto const & arc : mularcs) { 
            if (processed.count(arc.first) == 0 && graph->is_vec_T_consistent_color(arc.second)) {
              how_many_edges[arc.second].push_back(std::make_pair(vertex, arc.first));
            }             
          }
          processed.insert(vertex);
        }

        for (auto const & edges_with_color : how_many_edges) {
          if (edges_with_color.second.size() == 2) {
            auto const & first_edge = *edges_with_color.second.cbegin();
            auto const & second_edge = *(edges_with_color.second.cbegin() + 1);

            mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(first_edge.first);
            mularcs_t mularcs_y = graph->get_adjacent_multiedges_with_info(first_edge.second);

            twobreak_t twobreak;
            //KILL ME FOR THIS CODE, but I write it very fast and with end of working day. 
            if (mularcs_x.defined(second_edge.first) && !mularcs_x.defined(second_edge.second)
              && !mularcs_y.defined(second_edge.first) && !mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge, second_edge, edges_with_color.first);
            } else if (!mularcs_x.defined(second_edge.first) && mularcs_x.defined(second_edge.second)
              && !mularcs_y.defined(second_edge.first) && !mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge.first, first_edge.second, second_edge.second, second_edge.first, edges_with_color.first);
            } else if (!mularcs_x.defined(second_edge.first) && !mularcs_x.defined(second_edge.second)
              && mularcs_y.defined(second_edge.first) && !mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge.second, first_edge.first, second_edge.first, second_edge.second, edges_with_color.first);
            } else if (!mularcs_x.defined(second_edge.first) && !mularcs_x.defined(second_edge.second)
              && !mularcs_y.defined(second_edge.first) && mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge.second, first_edge.first, second_edge.second, second_edge.first, edges_with_color.first);
            } else if (mularcs_x.defined(second_edge.first) && !mularcs_x.defined(second_edge.second)
              && !mularcs_y.defined(second_edge.first) && mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge, second_edge, edges_with_color.first);
            } else if (!mularcs_x.defined(second_edge.first) && mularcs_x.defined(second_edge.second)
              && mularcs_y.defined(second_edge.first) && !mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge.first, first_edge.second, second_edge.second, second_edge.first, edges_with_color.first);
            }

            //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
            graph->apply(twobreak);
            ++number_rear;
            break;
          }
        }
      }

      if (number_rear != 0) { 
      isChanged = true;
    }
  } while (number_rear > 0); 

  return isChanged;   
}

/*
template<class graph_t>
bool Algorithm<graph_t>::stage5_3() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0;
      // go over all T-consistent multicolors
    for(auto ic = graph->cbegin_vec_T_consistent_color(); ic != graph->cend_vec_T_consistent_color(); ++ic) {
      const auto& Q = *ic;
      bool repeat = true;
      while(repeat) {
        repeat = false;
        std::map<vertex_t, std::set<arc_t> > EC; // reg. edges between diff. connected components of color Q
        std::map<vertex_t, std::set<arc_t> > EI; // irreg. edges of color Q
        utility::equivalence<vertex_t> CC = split_on_components(EC, EI, Q); 
    
        std::unordered_set<vertex_t> processed;
        // reg. edges between diff. connected components of color QQ
        for (const auto &reg_edge : EC) {
          if (processed.count(reg_edge.first) != 0) { 
            continue;
          }

          if (EI[reg_edge.first].size() != 0 || (reg_edge.second.size() % 2 != 0)) {
            continue;
          }

          std::set<vertex_t> verteces;
          partgraph_t bridges;
          for (auto const & arcs : reg_edge.second) { 
            verteces.insert(arcs.second);
            bridges.insert(arcs.first, arcs.second);
          }

          bool found = true;
          std::unordered_set<vertex_t> marks;
          std::list<twobreak_t> breaks;
          std::cerr << "NEW COMPONENT" << std::endl;
          for (auto const & arcs : reg_edge.second) {
            if (marks.count(arcs.second) == 0) {
              std::cerr << "Start to check " << arcs.second << std::endl;
              marks.insert(arcs.second);

              mularcs_t mularcs = graph->get_adjacent_multiedges(arcs.second);
              
              vertex_t match = ""; 
              size_t count_match = 0;
              for (auto marcs = mularcs.cbegin(); marcs != mularcs.cend(); ++marcs) {
                if (verteces.count(marcs->first) != 0) { 
                  match = marcs->first;
                  ++count_match;
                }
              }

              if (count_match == 1) {
                count_match = 0;
                mularcs_t mularcs_v = graph->get_adjacent_multiedges(match);
                for (auto marcs = mularcs_v.cbegin(); marcs != mularcs_v.cend(); ++marcs) {
                  if (verteces.count(marcs->first) != 0) { 
                    ++count_match;
                  }
                }       

                if (count_match == 1) {
                  marks.insert(match);
                  assert(bridges.find(match) != bridges.end());
                  vertex_t u = bridges.find(match)->second;
                  std::cerr << "perfect find " << u << " " << match << std::endl;
                  twobreak_t br2(arcs.first, arcs.second, u, match, Q);
                  breaks.push_back(br2);
                } else { 
                  found = false;
                }
              } else {
                found = false;
              }
            }

            if (!found) {
              break;
            }
          }

          if (found) {
            std::cerr << "ALL IS GOOD, do 2break" << std::endl;
            for (auto const & br2 : breaks) {
              graph->apply(br2);
              ++number_rear;
            }
            repeat = true;
            break;
          }
          std::cerr << "END COMPONENT" << std::endl;
        }
      }
    }

    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 

  return isChanged;   
}
*/

/*
template<class graph_t>
bool Algorithm<graph_t>::stage71() {
  size_t number_rear = 0; // number of rearrangements 

  for (auto const & x: *graph) {  
    mularcs_t const & mularcs = graph->get_adjacent_multiedges(x);
      
    bool found = false;
    for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
      vertex_t const & y = im->first; // Q == im->second - color of central edge
    
      if (y == Infty) {
        continue;
      }

      if (!graph->is_T_consistent_color(im->second)) { 
      mularcs_t mularcs_x = graph->get_adjacent_multiedges(x);
      mularcs_x.erase(y);

      mularcs_t mularcs_y = graph->get_adjacent_multiedges(y);
      mularcs_y.erase(x);
 
      size_t count = 0;
      bool is_y = false;
      std::pair<vertex_t, mcolor_t> pr;
      for (auto left = mularcs_x.begin(); left != mularcs_x.end(); ++left) {  
        if (!graph->is_T_consistent_color(left->second)) {
          pr = *left;
          ++count;
        } 
      }     

      for (auto right = mularcs_y.begin(); right != mularcs_y.end(); ++right) { 
        if (!graph->is_T_consistent_color(right->second)) {
          pr = *right;
          is_y = true;
          ++count;
        } 
      }

      if (count == 1 && pr.first != Infty) { 
        size_t count_f = graph->split_color(im->second).size(); 
        size_t count_s = graph->split_color(pr.second).size();

        if (count_f == count_s) { 
    //std::cerr << "Situation " << x << " " << y << " " << pr.first << std::endl;
          mularcs_t mul_x = graph->get_adjacent_multiedges_with_info(x, false);
          mularcs_t mul_y = graph->get_adjacent_multiedges_with_info(y, false);
          size_t count_central = calculate_cost(y, mul_x, mul_y);
          size_t count_another = 0;
          if (is_y) { 
            auto temp = graph->get_adjacent_multiedges_with_info(pr.first, false);
            temp.erase(y);
            count_another = calculate_cost(pr.first, graph->get_adjacent_multiedges_with_info(y, false), temp);
    } else {
            auto temp = graph->get_adjacent_multiedges_with_info(pr.first, false);
            temp.erase(y);
           count_another = calculate_cost(pr.first, graph->get_adjacent_multiedges_with_info(x, false), temp);
          } 

          if (count_central == std::min(count_central, count_another)) {
            mul_x.erase(y);
            std::vector<twobreak_t> history;  
            bool good = true;
            for (auto const &arc : mul_x) {
        vertex_t const & v = mularcs_y.get_vertex(arc.second); 
              if (!v.empty() && graph->is_vec_T_consistent_color(arc.second)) {
                  //std::cerr << "Two_break " << x << " " << arc.first << " " << y << " " << v << " " << genome_match::mcolor_to_name(arc.second) << std::endl;                
                  history.push_back(twobreak_t(x, arc.first, y, v, arc.second));
            //graph->apply(twobreak_t(x, arc.first, y, v, arc.second));
              } else { 
                good = false;
              }  
            }

            if (good) { 
              found = true;
              for (auto const & break2 : history) {
                graph->apply(break2);
                ++number_rear; 
              }
            } 
          } 
        }
      }
      }  
    } 
  }
 
  return (number_rear != 0);
}
*/

#endif