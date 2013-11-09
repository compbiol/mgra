#ifndef MBGRAPH_COLOR_H_
#define MBGRAPH_COLOR_H_

#include "mbgraph.h"

#include "genome_match.h"

template<class mcolor_t>
struct mbgraph_with_colors: public MBGraph { 
  typedef structure::Mularcs<mcolor_t> mularcs_t;
  typedef typename std::set<mcolor_t>::const_iterator citer; 
  
  template<class conf_t>
  mbgraph_with_colors(const std::vector<MBGraph::genome_t>& genomes, const conf_t& cfg); 

  mularcs_t get_adjacent_multiedges_with_info(const vertex_t& u, bool with_bad_edge = true);
  
  inline mcolor_t get_complement_color(const mcolor_t& color) { 
    assert(color.is_one_to_one_match()); 
    if (compliment_colors.count(color) == 0) { 
      const mcolor_t& temp = compute_complement_color(color);  
      compliment_colors.insert(std::make_pair(color, temp));
      compliment_colors.insert(std::make_pair(temp, color));
      return temp;
    }  
    return compliment_colors.find(color)->second;
  } 
  
  inline mcolor_t get_min_complement_color(const mcolor_t& color) {
    const mcolor_t& temp = get_complement_color(color);
    if (temp.size() > color.size() || (temp.size() == color.size() && temp > color)) {	
      return color;
    } 
    return temp;
  }	

  inline void registrate_viewed_edge(const vertex_t &u, const vertex_t& v) {
    not_mobile_edges.insert(u, v);
  } 

  inline void update_number_of_splits(size_t ns) {
    if (ns == 3) { 
      number_of_splits = this->count_local_graphs() + 1; 
    } else { 
      number_of_splits = ns;  
    } 
  } 

  mularcs_t get_adjacent_multiedges(const vertex_t& u) const; 
  size_t max_degree_split_color(const mcolor_t& color) const;

  bool is_simple_vertex(const vertex_t& v) const;
  bool is_indel_vertex(const vertex_t& v) const;  
  bool is_duplication_vertex(const vertex_t& v) const;
  bool is_have_self_loop(const vertex_t& v) const;

  std::map<vertex_t, std::set<vertex_t> > split_on_components(bool not_drop_complete_edge = true) const;
  bool are_adjacent_branches(const mcolor_t& A, const mcolor_t & B) const;

  inline mcolor_t get_complete_color() const {
    return complete_color;
  }

  inline bool is_T_consistent_color(const mcolor_t& color) const { 
    return (T_consistent_colors.count(color) > 0);
  } 

  inline size_t count_vec_T_consitent_color() const { 
    return vec_T_consistent_colors.size();
  } 
  
  inline bool is_vec_T_consistent_color(const mcolor_t& color) const {
    return (vec_T_consistent_colors.count(color) > 0);
  }

  inline citer cbegin_Tconsistent_color() const { //FIXME change
    return T_consistent_colors.cbegin(); 
  } 

  inline citer cend_Tconsistent_color() const { //FIXME change 
    return T_consistent_colors.cend(); 
  } 

  inline citer cbegin_T_consistent_color() const { //FIXME change
    return vec_T_consistent_colors.cbegin(); 
  } 

  inline citer cend_T_consistent_color() const { //FIXME change 
    return vec_T_consistent_colors.cend(); 
  } 

  inline mcolor_t get_root_color() const {
    return remove_color;
  } 
private: 
  inline mcolor_t compute_complement_color(const mcolor_t& color) const {
    mcolor_t answer; 
    for(size_t j = 0; j < count_local_graphs(); ++j) { 
      if (!color.defined(j)) { 
	answer.insert(j);
      } 
    } 
    return answer;
  }

  std::set<mcolor_t> split_color(const mcolor_t& color);
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

template<class mcolor_t>
template<class conf_t>
mbgraph_with_colors<mcolor_t>::mbgraph_with_colors(const std::vector<MBGraph::genome_t>& genomes, const conf_t& cfg) 
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
    it->build_vec_T_consistent_colors(nodes_color);
  }

  for (const auto& color : nodes_color) {    
    if (color.includes(remove_color)) {
      vec_T_consistent_colors.insert(get_complement_color(color));
    } else {
      vec_T_consistent_colors.insert(color);
    } 
  } 

  vec_T_consistent_colors.erase(remove_color);
#endif

#ifdef VERSION2  
  std::set<mcolor_t> nodes_color;
  for (auto it = cfg.cbegin_trees(); it != cfg.cend_trees(); ++it) {
    it->build_vec_T_consistent_colors(nodes_color);
  }

  nodes_color.erase(complete_color);

  mcolor_t root_color = complete_color;
  size_t est = (complete_color.size() / 2 + complete_color.size() % 2);
 
  for (const auto& color : nodes_color) {
    if (color.size() >= est && root_color.size() >= color.size()) { 
      root_color = color;
    }   
  } 

  remove_color = root_color;
  
  for (const auto& color : nodes_color) {
    const auto& temp = compute_complement_color(color);   
    if (!color.includes(root_color)) {
      vec_T_consistent_colors.insert(color);
    } else if (!temp.includes(root_color)) { 
      vec_T_consistent_colors.insert(temp);
    } 
    T_consistent_colors.insert(color);
    T_consistent_colors.insert(temp);
  } 
#else 
  for (auto it = cfg.cbegin_trees(); it != cfg.cend_trees(); ++it) {
    it->build_vec_T_consistent_colors(vec_T_consistent_colors);
  }
  vec_T_consistent_colors.erase(complete_color);
#endif

  vec_T_consistent_colors.erase(cfg.get_target());
  
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
   
  for (const auto &vtc : vec_T_consistent_colors) {
    T_consistent_colors.insert(vtc);
    T_consistent_colors.insert(compute_complement_color(vtc));
  }
}

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::is_simple_vertex(const vertex_t& v) const {
  const mularcs_t& mularcs = get_adjacent_multiedges(v);
  if (mularcs.size() == 2 && mularcs.cbegin()->second.is_one_to_one_match() && mularcs.crbegin()->second.is_one_to_one_match() 
      && !is_duplication_vertex(v) && !is_indel_vertex(v)) {
	return true;  
  }
  return false; 
}  

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::is_have_self_loop(const vertex_t& v) const {
  const mularcs_t& mularcs = get_adjacent_multiedges(v);
  return mularcs.defined(v);
} 
 
template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::is_indel_vertex(const vertex_t& v) const {
  if (is_duplication_vertex(v)) {
    return false; 
  }  
 
  const mcolor_t& un = get_adjacent_multiedges(v).union_multicolors();
	
  if (!un.is_one_to_one_match() || (un == complete_color)) {
    return false; 
  }  

  return true; 
}

template<class mcolor_t>  
bool mbgraph_with_colors<mcolor_t>::is_duplication_vertex(const vertex_t& v) const {	
  if (is_have_self_loop(v)) {
    return true;
  } 

  const mularcs_t& mularcs = get_adjacent_multiedges(v);
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

template<class mcolor_t>
structure::Mularcs<mcolor_t> mbgraph_with_colors<mcolor_t>::get_adjacent_multiedges(const vertex_t& u) const { 
  if (u == Infty) {
    std::cerr << "mularcs ERROR: Infinite input" << std::endl;
    exit(1);
  }

  mularcs_t output;
  for (size_t i = 0; i < count_local_graphs(); ++i) {
    auto iters = local_graph[i].equal_range(u);
    for (auto it = iters.first; it != iters.second; ++it) { 
      output.insert(it->second, i); 
    }
  }

  return output;
} 

template<class mcolor_t>
structure::Mularcs<mcolor_t> mbgraph_with_colors<mcolor_t>::get_adjacent_multiedges_with_info(const vertex_t& u, bool with_bad_edge) { 
  mularcs_t&& output = get_adjacent_multiedges(u);
  
  if (number_of_splits != 1) { 
    mularcs_t split; 
    for(const auto &arc : output) {
      if ((!with_bad_edge || (with_bad_edge && !not_mobile_edges.defined(u, arc.first))) 
	   && !is_vec_T_consistent_color(arc.second) && arc.second.size() < count_local_graphs()) {
	const auto& colors = split_color(arc.second);
	for(const auto &color : colors) {
	  split.insert(arc.first, color); 
	}
      } else { 
	split.insert(arc.first, arc.second); 
      }
    }
    return split; 
  }

  return std::move(output);
} 

template<class mcolor_t>
size_t mbgraph_with_colors<mcolor_t>::max_degree_split_color(const mcolor_t& color) const {
  if (is_vec_T_consistent_color(color)) {
    return 1; 
  } else { 
    utility::equivalence<size_t> equiv;
    std::for_each(color.cbegin(), color.cend(), [&] (const std::pair<size_t, size_t>& col) -> void {
      equiv.addrel(col.first, col.first);
    }); 

    for (const auto &vtc: vec_T_consistent_colors) { 
      mcolor_t inter_color(vtc, color, mcolor_t::Intersection);
      if (inter_color.size() >= 2 && inter_color.size() == vtc.size()) {
        std::for_each(inter_color.cbegin(), inter_color.cend(), [&] (const std::pair<size_t, size_t>& col) -> void {
          equiv.addrel(col.first, inter_color.cbegin()->first);
        });
      }
    }

    equiv.update();

    return equiv.classes();
  } 
}

template<class mcolor_t>
std::set<mcolor_t> mbgraph_with_colors<mcolor_t>::split_color(const mcolor_t& color) {
  if (is_vec_T_consistent_color(color) || (number_of_splits == 1)) {
    return std::set<mcolor_t>({color}); 
  } else {  
    std::set<mcolor_t> answer;
    
    /*if (hashing_split_colors.count(std::make_pair(color, number_of_splits)) != 0) {
      answer = hashing_split_colors.find(std::make_pair(color, number_of_splits))->second;
    } else {   */
      utility::equivalence<size_t> equiv;
      std::for_each(color.cbegin(), color.cend(), [&] (const std::pair<size_t, size_t>& col) -> void {
        equiv.addrel(col.first, col.first);
      }); 

      bool tcons = false; 
      #ifdef VERSION2  
      for (const auto &tc: T_consistent_colors) { 
        mcolor_t inter_color(tc, color, mcolor_t::Intersection);
        if (inter_color.size() >= 2 && inter_color.size() == tc.size()) { 
  	  tcons = (vec_T_consistent_colors.count(color) == 0);
          std::for_each(inter_color.cbegin(), inter_color.cend(), [&] (const std::pair<size_t, size_t>& col) -> void {
            equiv.addrel(col.first, inter_color.cbegin()->first);
          });
        }
      }
      #else 
        for (const auto &vtc: vec_T_consistent_colors) { 
          mcolor_t inter_color(vtc, color, mcolor_t::Intersection);
          if (inter_color.size() >= 2 && inter_color.size() == vtc.size()) {
            std::for_each(inter_color.cbegin(), inter_color.cend(), [&] (const std::pair<size_t, size_t>& col) -> void {
              equiv.addrel(col.first, inter_color.cbegin()->first);
            });
          }
        }
      }
      #endif

      if (tcons) { 
	answer.insert(color);         
      } else { 
        equiv.update();
        const std::map<size_t, mcolor_t>& classes = equiv.get_eclasses<mcolor_t>(); 
        //std::cerr << "color " << genome_match::mcolor_to_name(color) << std::endl;
        for(const auto &col : classes) {
          answer.insert(col.second);
        }

        if (answer.size() > number_of_splits) { 
          answer.clear(); 
          answer.insert(color);  	
        }
      }  

      //hashing_split_colors.insert(std::make_pair(std::make_pair(color, number_of_splits), answer));
    //}
    return answer;
  }
} 

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::are_adjacent_branches(const mcolor_t& mcolor_a, const mcolor_t & mcolor_b) const {
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

template<class mcolor_t>
std::map<vertex_t, std::set<vertex_t> > mbgraph_with_colors<mcolor_t>::split_on_components(bool not_drop_complete_edge) const { 
  utility::equivalence<vertex_t> connected_components; // connected components

  for(const auto &x : vertex_set) {
    if (!not_drop_complete_edge) { 
	connected_components.addrel(x, x);
    } 

    const mularcs_t& mularcs = get_adjacent_multiedges(x); 

    if (not_drop_complete_edge && mularcs.number_unique_edge() == 1 && mularcs.union_multicolors() == get_complete_color()) { 
      continue; // ignore complete multiedges
    } 

    std::for_each(mularcs.cbegin(), mularcs.cend(), [&] (const std::pair<vertex_t, mcolor_t>& arc) -> void {    
      if (arc.first != Infty) { 
	connected_components.addrel(x, arc.first);
      } 
    });
  }
		    
  connected_components.update();
   
  return connected_components.get_eclasses<std::set<vertex_t> >(); 
}

#endif

