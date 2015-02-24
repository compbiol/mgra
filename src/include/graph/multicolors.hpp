#ifndef MULTICOLORS_HPP
#define MULTICOLORS_HPP

template<class mcolor_t>
struct Multicolors { 
  
  Multicolors(); 

  /*
   * Split input multiolor COLOR on T-consistent multicolors according number of splits.
   * If multicolor split to more than NUMBER_SPLITS value return input color. 
   * If less returned set of T-consistent color.
   */
  std::set<mcolor_t> split_color_on_tc_color(mcolor_t const & color, size_t number_splits) const;

  /*
   * Split input multiolor COLOR on vec{T}-consistent multicolors.
   */
  std::set<mcolor_t> split_color_on_vtc_color(mcolor_t const & color) const;
  std::set<mcolor_t> split_color_on_next_vtc_color(mcolor_t const & color) const;

  mcolor_t get_min_addit_color_for_tc(mcolor_t const & color) const;
  
  inline mcolor_t const & get_complement_color(mcolor_t const & color) { 
    assert(color.is_one_to_one_match()); 
    if (complement_colors.count(color) == 0) {  
      mcolor_t const & temp = compute_complement_color(color);  
      complement_colors.insert(std::make_pair(color, temp));
      complement_colors.insert(std::make_pair(temp, color));
    }  
    return complement_colors.find(color)->second;
  } 
  
  inline bool is_T_consistent_color(mcolor_t const & color) const { 
    return (T_consistent_colors.count(color) > 0);
  } 

  inline bool is_vec_T_consistent_color(mcolor_t const & color) const {
    return (vec_T_consistent_colors.count(color) > 0);
  }

  std::vector<std::tuple<mcolor_t, mcolor_t, mcolor_t> > get_medians_colors() const { 
    return median_colors;
  }

  DECLARE_GETTER(mcolor_t const &, complete_color, complete_color)
  DECLARE_GETTER(mcolor_t const &, remove_color, root_color)
  DECLARE_DELEGATE_CONST_METHOD( size_t, vec_T_consistent_colors, count_vec_T_consitent_color, size )

  using citer = typename std::set<mcolor_t>::const_iterator; 
  DECLARE_CONST_ITERATOR( citer, T_consistent_colors, cbegin_T_consistent_color, cbegin )  
  DECLARE_CONST_ITERATOR( citer, T_consistent_colors, cend_T_consistent_color, cend )
  DECLARE_CONST_ITERATOR( citer, vec_T_consistent_colors, cbegin_vec_T_consistent_color, cbegin )  
  DECLARE_CONST_ITERATOR( citer, vec_T_consistent_colors, cend_vec_T_consistent_color, cend )
  
private: 
  using tree_t = typename structure::BinaryTree<mcolor_t>; 
  using node_t = typename tree_t::Node; 
  using colors_median_t = std::tuple<mcolor_t, mcolor_t, mcolor_t>; 

  /*
   *
   */
  void get_vector_colors_from_tree(std::unique_ptr<node_t> const & current, 
                                          std::set<mcolor_t>& vec_colors) const; 
  void get_median_colors_from_tree(std::unique_ptr<node_t> const & current); 

  /*
   *
   */
  inline mcolor_t compute_complement_color(mcolor_t const & color) const {
    mcolor_t answer(complete_color, color, mcolor_t::Difference); 
    return answer;
  }

protected:
  mcolor_t complete_color;
  mcolor_t remove_color; 
  
  std::set<mcolor_t> T_consistent_colors;
  std::set<mcolor_t> vec_T_consistent_colors;
  
  std::map<mcolor_t, mcolor_t> complement_colors;
  
  std::vector<colors_median_t> median_colors;
}; 

template<class mcolor_t>
Multicolors<mcolor_t>::Multicolors() {
  //Leaf have a vec{T}-consistent multicolor.
  for (size_t i = 0; i < cfg::get().get_count_genomes(); ++i) {
    complete_color.insert(i);
    vec_T_consistent_colors.insert(mcolor_t(i));
  }

  // Get all vec{T}-consistent colors corresponding input [sub]tree[s]. 
  std::set<mcolor_t> nodes_color;
  for (auto const & tree : cfg::get().phylotrees) {
    get_vector_colors_from_tree(tree.get_root(), nodes_color); 
  }
  nodes_color.erase(complete_color);

  // Get all median colors corresponding input [sub]tree[s]. 
  for (auto const & tree : cfg::get().phylotrees) {
    get_median_colors_from_tree(tree.get_root());
  }
  
  //If target is empty we put root in nearest node. Work fine only complete tree.
  //Need tested on subtrees. 
  if (cfg::get().how_build != target_algo) { 
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
    } 
  } else { 
    remove_color = cfg::get().target_mcolor;
    nodes_color.erase(cfg::get().target_mcolor);
  }
  vec_T_consistent_colors = nodes_color;

  //check consistency for multicolors
  for (auto id = vec_T_consistent_colors.cbegin(); id != vec_T_consistent_colors.cend(); ++id) {
    for(auto jd = id; jd != vec_T_consistent_colors.end(); ++jd) {
      mcolor_t color(*id, *jd, mcolor_t::Intersection);
      if (!color.empty() && color.size() != id->size() && color.size() != jd->size()) {
      	vec_T_consistent_colors.erase(jd++);
      	--jd;
      }
    }
  }
   
  //init T consistent multicolors (vec{T}-consistent multicolors and addition).  
  for (auto const & vtc : vec_T_consistent_colors) {
    T_consistent_colors.insert(vtc);
    T_consistent_colors.insert(compute_complement_color(vtc));
  }
}

/*
 *
 */
template<class mcolor_t>
void Multicolors<mcolor_t>::get_vector_colors_from_tree(std::unique_ptr<node_t> const & current,
                                     std::set<mcolor_t>& vec_colors) const { 
  if (current->get_left_child()) {
    get_vector_colors_from_tree(current->get_left_child(), vec_colors); 
  }

  vec_colors.insert(current->get_data());

  if (current->get_right_child()) {
    get_vector_colors_from_tree(current->get_right_child(), vec_colors); 
  }
}

template<class mcolor_t>
void Multicolors<mcolor_t>::get_median_colors_from_tree(std::unique_ptr<node_t> const & current) {
  auto const & left = current->get_left_child();
  if (left) { 
    get_median_colors_from_tree(left); 
  }

  auto const & right = current->get_right_child();
  if (right) { 
    get_median_colors_from_tree(right); 
  }

  if (left && right && current->get_parent() != nullptr) { 
    mcolor_t parent = mcolor_t(complete_color, mcolor_t(left->get_data(), right->get_data(), mcolor_t::Union), mcolor_t::Difference);
    median_colors.push_back(std::make_tuple(left->get_data(), right->get_data(), parent));
  }
}

/*
 *
 */
template<class mcolor_t>
std::set<mcolor_t> Multicolors<mcolor_t>::split_color_on_tc_color(mcolor_t const & color, size_t number_splits) const {
  std::set<mcolor_t> answer;    
  
  if (is_T_consistent_color(color) || (number_splits == 1)) {
    answer.insert(color);
  } else {  
    
    /*if (hashing_split_colors.count(std::make_pair(color, number_of_splits)) != 0) {
      answer = hashing_split_colors.find(std::make_pair(color, number_of_splits))->second;
    } else {   */
    utility::equivalence<size_t> equiv;
    std::for_each(color.cbegin(), color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
      equiv.addrel(col.first, col.first);
    }); 

    for (auto const & tc: T_consistent_colors) { 
      mcolor_t inter_color(tc, color, mcolor_t::Intersection);
      if (inter_color.size() >= 2 && inter_color.size() == tc.size()) { 
        std::for_each(inter_color.cbegin(), inter_color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
          equiv.addrel(col.first, inter_color.cbegin()->first);
        });
      }
    }
    
    equiv.update();
    std::map<size_t, mcolor_t> const & classes = equiv.get_eclasses<mcolor_t>(); 
    //std::cerr << "color " << genome_match::mcolor_to_name(color) << std::endl;
    for(auto const & col : classes) {
      answer.insert(col.second);
    }

    if (answer.size() > number_splits) { 
      answer.clear(); 
      answer.insert(color);   
    }
  }  
  return answer;  
}

template<class mcolor_t>
std::set<mcolor_t> Multicolors<mcolor_t>::split_color_on_vtc_color(mcolor_t const & color) const {
  std::set<mcolor_t> answer; 

  if (is_vec_T_consistent_color(color)) {
    answer.insert(color);
  } else {  
    std::set<mcolor_t> answer;    

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
    std::map<size_t, mcolor_t> const & classes = equiv.get_eclasses<mcolor_t>(); 
    for(auto const & col : classes) {
      answer.insert(col.second);
    }
  }

  return answer;
}

template<class mcolor_t>
std::set<mcolor_t> Multicolors<mcolor_t>::split_color_on_next_vtc_color(mcolor_t const & color) const { 
  assert(is_vec_T_consistent_color(color));
  std::set<mcolor_t> answer;    

  utility::equivalence<size_t> equiv;
  std::for_each(color.cbegin(), color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
    equiv.addrel(col.first, col.first);
  }); 
  
  for (auto const & vtc: vec_T_consistent_colors) { 
    if (vtc != color) {  
      mcolor_t inter_color(vtc, color, mcolor_t::Intersection);
      if (inter_color.size() >= 2 && inter_color.size() == vtc.size()) {
        std::for_each(inter_color.cbegin(), inter_color.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
          equiv.addrel(col.first, inter_color.cbegin()->first);
        });
      }
    } 
  }

  equiv.update();
  std::map<size_t, mcolor_t> const & classes = equiv.get_eclasses<mcolor_t>(); 
  for(auto const & col : classes) {
    answer.insert(col.second);
  }

  if (answer.size() == 1 && *answer.begin() == color) { 
    answer.clear();
  } 
  
  return answer;
}

template<class mcolor_t>
mcolor_t Multicolors<mcolor_t>::get_min_addit_color_for_tc(mcolor_t const & color) const { 
  mcolor_t min_color = get_complete_color();
  for (auto col = cbegin_T_consistent_color(); col != cend_T_consistent_color(); ++col) {
    if (col->includes(color)) {
      mcolor_t diff_color(*col, color, mcolor_t::Difference);
      if (diff_color.size() < min_color.size()) {
        min_color = diff_color;
      } 
    } 
  } 
  return min_color;
}


#endif

