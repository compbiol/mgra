#ifndef MBGRAPH_COLOR_H_
#define MBGRAPH_COLOR_H_

#include "graph/mbgraph.h"

#include "genome_match.h"

template<class mcolor_type>
struct mbgraph_with_colors: public MBGraph { 
  typedef mcolor_type mcolor_t;
  typedef structure::Mularcs<mcolor_t> mularcs_t;
  typedef typename std::set<mcolor_t>::const_iterator citer; 
  
  template<class conf_t>
  mbgraph_with_colors(std::vector<MBGraph::genome_t> const & genomes, conf_t const & cfg); 

  mularcs_t get_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge = true);
  
  inline mcolor_t get_complement_color(mcolor_t const & color) { 
    assert(color.is_one_to_one_match()); 
    if (compliment_colors.count(color) == 0) {  
      mcolor_t const & temp = compute_complement_color(color);  
      compliment_colors.insert(std::make_pair(color, temp));
      compliment_colors.insert(std::make_pair(temp, color));
      return temp;
    }  
    return compliment_colors.find(color)->second;
  } 
  
  inline mcolor_t get_min_complement_color(mcolor_t const & color) {
    mcolor_t const & temp = get_complement_color(color);
    if (temp.size() > color.size() || (temp.size() == color.size() && temp > color)) {	
      return color;
    } 
    return temp;
  }	

  inline void registrate_viewed_edge(vertex_t const & u, vertex_t const & v) {
    not_mobile_edges.insert(u, v);
  } 

  inline void update_number_of_splits(size_t ns) {
    if (ns == 3) { 
      number_of_splits = this->count_local_graphs() + 1; 
    } else { 
      number_of_splits = ns;  
    } 
  } 

  mularcs_t get_adjacent_multiedges(vertex_t const & u) const; 
  size_t max_degree_split_color(mcolor_t const & color) const;

  bool is_simple_vertex(vertex_t const & v) const;
  bool is_indel_vertex(vertex_t const & v) const;  
  bool is_duplication_vertex(vertex_t const & v) const;
  bool is_have_self_loop(vertex_t const & v) const;

  std::map<vertex_t, std::set<vertex_t> > split_on_components(bool not_drop_complete_edge = true) const;
  bool are_adjacent_branches(mcolor_t const & A, mcolor_t const & B) const;

  inline mcolor_t get_complete_color() const {
    return complete_color;
  }

  inline bool is_T_consistent_color(mcolor_t const & color) const { 
    return (T_consistent_colors.count(color) > 0);
  } 

  inline bool is_vec_T_consistent_color(mcolor_t const & color) const {
    return (vec_T_consistent_colors.count(color) > 0);
  }

  inline size_t count_vec_T_consitent_color() const { 
    return vec_T_consistent_colors.size();
  } 
  
  inline citer cbegin_Tconsistent_color() const { //FIXME change
    return T_consistent_colors.cbegin(); 
  } 

  inline citer cend_Tconsistent_color() const { //FIXME change 
    return T_consistent_colors.cend(); 
  } 

  inline citer cbegin_vec_T_consistent_color() const { 
    return vec_T_consistent_colors.cbegin(); 
  } 

  inline citer cend_vec_T_consistent_color() const { 
    return vec_T_consistent_colors.cend(); 
  } 

  inline mcolor_t const & get_root_color() const {
    return remove_color;
  } 
 
  std::set<mcolor_t> split_color(mcolor_t const & color);

private: 
  inline mcolor_t compute_complement_color(mcolor_t const & color) const {
    mcolor_t answer; 
    for(size_t j = 0; j < count_local_graphs(); ++j) { 
      if (!color.defined(j)) { 
        answer.insert(j);
      } 
    } 
    return answer;
  }

protected:
  size_t number_of_splits; 
  mcolor_t complete_color;
  mcolor_t remove_color; 
  std::map<mcolor_t, mcolor_t> compliment_colors;
  std::set<mcolor_t> T_consistent_colors;
  std::set<mcolor_t> vec_T_consistent_colors;
  edges_t not_mobile_edges;
  std::map<std::pair<mcolor_t, size_t>, std::set<mcolor_t> > hashing_split_colors;
}; 

template<class mcolor_type>
template<class conf_t>
mbgraph_with_colors<mcolor_type>::mbgraph_with_colors(std::vector<MBGraph::genome_t> const & genomes, conf_t const & cfg) 
: MBGraph(genomes)
, number_of_splits(1)
{
  for (size_t i = 0; i < genomes.size(); ++i) {
    complete_color.insert(i);
    vec_T_consistent_colors.insert(mcolor_t(i));
  }

#ifdef ROOT_LEAF
  remove_color = mcolor_t(genomes.size() - 1);
  
  std::set<mcolor_t> nodes_color;
  for (auto it = cfg.cbegin_trees(); it != cfg.cend_trees(); ++it) {
    auto const & tree_color = it->build_vec_T_consistent_colors();
    nodes_color.insert(tree_color.begin(), tree_color.end());
  }

  for (auto const & color : nodes_color) {    
    if (color.includes(remove_color)) {
      vec_T_consistent_colors.insert(get_complement_color(color));
    } else {
      vec_T_consistent_colors.insert(color);
    } 
  } 

  vec_T_consistent_colors.erase(remove_color);
#endif

  std::set<mcolor_t> nodes_color;
  for (auto it = cfg.cbegin_trees(); it != cfg.cend_trees(); ++it) {
    auto const & tree_color = it->build_vec_T_consistent_colors();
    nodes_color.insert(tree_color.cbegin(), tree_color.cend());
  }
  nodes_color.erase(complete_color);

#ifndef VERSION1  
  if (cfg.get_target().empty()) { 
    //mcolor_t root_color = complete_color;
    //size_t est = (complete_color.size() / 2 + complete_color.size() % 2);
 
    for (auto const & color : nodes_color) {
      auto const & compl_color = compute_complement_color(color);
      if (nodes_color.find(compl_color) != nodes_color.end()) {
        if (color.size() > compl_color.size()) { 
          remove_color = color;
        } else { 
          remove_color = compl_color; 
        } 
        nodes_color.erase(remove_color);
        break; 
      } 
      //if (color.size() >= est && root_color.size() >= color.size()) { 
        //root_color = color;
      //}   
    } 

    //remove_color = root_color;
    //nodes_color.erase(remove_color);
    //nodes_color.insert(compute_complement_color(remove_color));
  }
#endif 

  vec_T_consistent_colors = nodes_color;
  if (!cfg.get_target().empty()) { 
    remove_color = cfg.get_target();
    //vec_T_consistent_colors.insert(remove_color);
    vec_T_consistent_colors.erase(cfg.get_target());
  } 

  //check consistency
  for (auto id = vec_T_consistent_colors.cbegin(); id != vec_T_consistent_colors.cend(); ++id) {
    for(auto jd = id; jd != vec_T_consistent_colors.end(); ++jd) {
      mcolor_t color(*id, *jd, mcolor_t::Intersection);
      if (!color.empty() && color.size() != id->size() && color.size() != jd->size()) {
      	vec_T_consistent_colors.erase(jd++);
      	--jd;
      }
    }
  }
   
  for (auto const & vtc : vec_T_consistent_colors) {
    T_consistent_colors.insert(vtc);
    T_consistent_colors.insert(compute_complement_color(vtc));
  }
}

template<class mcolor_type>
bool mbgraph_with_colors<mcolor_type>::is_simple_vertex(vertex_t const & v) const {
  mularcs_t const & mularcs = get_adjacent_multiedges(v);
  if (mularcs.size() == 2 && mularcs.cbegin()->second.is_one_to_one_match() && mularcs.crbegin()->second.is_one_to_one_match() && !is_duplication_vertex(v) && !is_indel_vertex(v)) {
	return true;  
  }
  return false; 
}  

template<class mcolor_type>
bool mbgraph_with_colors<mcolor_type>::is_have_self_loop(vertex_t const & v) const {
  mularcs_t const & mularcs = get_adjacent_multiedges(v);
  return mularcs.defined(v);
} 
 
template<class mcolor_type>
bool mbgraph_with_colors<mcolor_type>::is_indel_vertex(vertex_t const & v) const {
  if (is_duplication_vertex(v)) {
    return false; 
  }  
 
  mcolor_t const & un = get_adjacent_multiedges(v).union_multicolors();
	
  if (!un.is_one_to_one_match() || (un == complete_color)) { //FIXME
    return false; 
  }  

  return true; 
}

template<class mcolor_type>  
bool mbgraph_with_colors<mcolor_type>::is_duplication_vertex(vertex_t const & v) const {	
  if (is_have_self_loop(v)) {
    return true;
  } 

  mularcs_t const & mularcs = get_adjacent_multiedges(v);
  for(auto im = mularcs.cbegin(); im != mularcs.cend(); ++im) { 
    if (!im->second.is_one_to_one_match()) {
	return true; 
    }

    for(auto it = mularcs.cbegin(); it != mularcs.cend(); ++it) {
      if (*im != *it) { 
        mcolor_t color(im->second, it->second, mcolor_t::Intersection);
        if (!color.empty()) { 
	  return true; 
        }
      }
    } 
  }  
  return false; 
} 

template<class mcolor_type>
structure::Mularcs<mcolor_type> mbgraph_with_colors<mcolor_type>::get_adjacent_multiedges(vertex_t const & u) const { 
  if (u == Infty) {
    std::cerr << "mularcs ERROR: Infinite input" << std::endl;
    exit(1);
  }

  mularcs_t output;
  for (size_t i = 0; i < count_local_graphs(); ++i) {
    auto iters = m_local_graphs[i].equal_range(u);
    for (auto it = iters.first; it != iters.second; ++it) { 
      output.insert(it->second, i); 
    }
  }

  return output;
} 

template<class mcolor_type>
structure::Mularcs<mcolor_type> mbgraph_with_colors<mcolor_type>::get_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge) { 
  mularcs_t output = get_adjacent_multiedges(u);
  
  if (number_of_splits != 1) { 
    mularcs_t split; 
    for(auto const & arc : output) {
      if ((!with_bad_edge || (with_bad_edge && !not_mobile_edges.defined(u, arc.first))) 
	   && !is_T_consistent_color(arc.second) && arc.second.size() < count_local_graphs()) {
        auto const & colors = split_color(arc.second);
        for(auto const & color : colors) {  
          split.insert(arc.first, color); 
        }
      } else { 
        split.insert(arc.first, arc.second); 
      }
    }
    return split; 
  }

  return output;
} 

template<class mcolor_type>
size_t mbgraph_with_colors<mcolor_type>::max_degree_split_color(mcolor_t const & color) const {
  if (is_vec_T_consistent_color(color)) {
    return 1; 
  } else { 
    utility::equivalence<size_t> equiv;
    std::for_each(color.cbegin(), color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
      equiv.addrel(col.first, col.first);
    }); 

    for (auto const & vtc: vec_T_consistent_colors) { 
      mcolor_t inter_color(vtc, color, mcolor_t::Intersection);
      if (inter_color.size() >= 2 && inter_color.size() == vtc.size()) {
        std::for_each(inter_color.cbegin(), inter_color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
          equiv.addrel(col.first, inter_color.cbegin()->first);
        });
      }
    }

    equiv.update();

    return equiv.classes();
  } 
}

template<class mcolor_type>
std::set<mcolor_type> mbgraph_with_colors<mcolor_type>::split_color(mcolor_t const & color) {
  if (is_T_consistent_color(color) || (number_of_splits == 1)) {
    return std::set<mcolor_t>({color}); 
  } else {  
    std::set<mcolor_t> answer;
    
    /*if (hashing_split_colors.count(std::make_pair(color, number_of_splits)) != 0) {
      answer = hashing_split_colors.find(std::make_pair(color, number_of_splits))->second;
    } else {   */
      utility::equivalence<size_t> equiv;
      std::for_each(color.cbegin(), color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
        equiv.addrel(col.first, col.first);
      }); 

      #ifndef VERSION1  
      for (auto const & tc: T_consistent_colors) { 
        mcolor_t inter_color(tc, color, mcolor_t::Intersection);
        if (inter_color.size() >= 2 && inter_color.size() == tc.size()) { 
          std::for_each(inter_color.cbegin(), inter_color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
            equiv.addrel(col.first, inter_color.cbegin()->first);
          });
        }
      }
      #else 
        for (auto const & vtc: vec_T_consistent_colors) { 
          mcolor_t inter_color(vtc, color, mcolor_t::Intersection);
          if (inter_color.size() >= 2 && inter_color.size() == vtc.size()) {
            std::for_each(inter_color.cbegin(), inter_color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
              equiv.addrel(col.first, inter_color.cbegin()->first);
            });
          }
        }
      }
      #endif

      equiv.update();
      std::map<size_t, mcolor_t> const & classes = equiv.get_eclasses<mcolor_t>(); 
      //std::cerr << "color " << genome_match::mcolor_to_name(color) << std::endl;
      for(auto const & col : classes) {
        answer.insert(col.second);
      }

      if (answer.size() > number_of_splits) { 
        answer.clear(); 
        answer.insert(color);  	
      }

      //hashing_split_colors.insert(std::make_pair(std::make_pair(color, number_of_splits), answer));
    //}
    return answer;
  }
} 

template<class mcolor_type>
bool mbgraph_with_colors<mcolor_type>::are_adjacent_branches(mcolor_t const & mcolor_a, mcolor_t const & mcolor_b) const {
  if (!is_T_consistent_color(mcolor_a) || !is_T_consistent_color(mcolor_b)) { 
    return false;
  } 

  mcolor_t color1; 
  mcolor_t color2;
  
  if (mcolor_a.size() >= mcolor_b.size()) {
    color1 = mcolor_a;
    color2 = mcolor_b;
  } else {
    color1 = mcolor_b;
    color2 = mcolor_a;
  }
  
  mcolor_t diff_color(color1, color2, mcolor_t::Difference);
  
  if (diff_color.size() == color1.size() - color2.size() && is_T_consistent_color(diff_color)) { 		
    return true;
  } 
  
  mcolor_t union_color(color1, color2, mcolor_t::Union); 
  if (union_color.size() == color1.size() + color2.size() && is_T_consistent_color(union_color)) { 	
    return true;
  } 
  
  return false;
}

template<class mcolor_type>
std::map<vertex_t, std::set<vertex_t> > mbgraph_with_colors<mcolor_type>::split_on_components(bool not_drop_complete_edge) const { 
  utility::equivalence<vertex_t> connected_components; // connected components

  for(auto const & x : vertex_set) {
    if (!not_drop_complete_edge) { 
      connected_components.addrel(x, x);
    } 

    mularcs_t const & mularcs = get_adjacent_multiedges(x); 

    if (not_drop_complete_edge && get_adjacent_multiedges(x).number_unique_edge() == 1 && mularcs.union_multicolors() == get_complete_color()) { 
      continue; // ignore complete multiedges
    } 

    std::for_each(mularcs.cbegin(), mularcs.cend(), [&] (std::pair<vertex_t, mcolor_t> const & arc) -> void {    
      if (arc.first != Infty) { 
        connected_components.addrel(x, arc.first);
      } 
    });
  }
		    
  connected_components.update();
   
  return connected_components.get_eclasses<std::set<vertex_t> >(); 
}

#endif

