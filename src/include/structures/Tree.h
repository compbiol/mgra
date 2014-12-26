#ifndef TREE_HPP 
#define TREE_HPP

namespace structure { 

template<class mcolor_t>
struct BinaryTree { 
  typedef std::set<mcolor_t> colors_median_t; 

  BinaryTree(std::string const & st, std::unordered_map<std::string, size_t> const & genome_number, std::vector<std::string> const & priority_name) 
  : root(new Node(nullptr, st, genome_number, priority_name))
  {
  } 

  void get_nodes(std::vector<std::string>& info) const {
    root->get_nodes(info);	
  }

  std::set<mcolor_t> build_vec_T_consistent_colors() const {
    std::set<mcolor_t> dicolor;
    root->get_dicolors(dicolor);
    return dicolor;	
  }

  std::map<mcolor_t, std::string> get_name_for_colors() const {
    std::map<mcolor_t, std::string> colors;
    root->get_name_for_colors(colors);
    return colors;	
  }

  std::vector<colors_median_t> get_median_colors() const { 
    std::vector<colors_median_t> medians; 
    root->get_median_colors(root->data, medians);
    return medians;
  }

  template<class linearizator_t>
  void walk_and_linearizeate(
    linearizator_t const & linearizator, 
    std::map<mcolor_t, typename linearizator_t::partgraph_t> & graphs, 
    std::map<mcolor_t, typename linearizator_t::transform_t> & transformations) const 
  {
    root->walk_and_linearizeate(linearizator, graphs, transformations);
  }

private: 
  struct Node { 
    std::string name;   
    std::vector<std::string> childs;
    
    mcolor_t data;

    Node * const parent;    
    std::unique_ptr<Node> left_child;    
    std::unique_ptr<Node> right_child; 

    Node() 
    : parent(nullptr) 
    , left_child (nullptr)
    , right_child (nullptr) 
    { 
    }

    Node(Node * const par, std::string const & tree, 
      std::unordered_map<std::string, size_t> const & genome_number, 
      std::vector<std::string> const & priority_name);     

    std::string get_nodes(std::vector<std::string>& info) const;

    void get_median_colors(mcolor_t const & complete_color, std::vector<colors_median_t>& medians) const; 
    void get_dicolors(std::set<mcolor_t>& dicolor) const; 
    void get_name_for_colors(std::map<mcolor_t, std::string>& colors) const;

    template<class linearizator_t>
    void walk_and_linearizeate(
      linearizator_t const & linearizator, 
      std::map<mcolor_t, typename linearizator_t::partgraph_t> & graphs, 
      std::map<mcolor_t, typename linearizator_t::transform_t> & transformations) const;
  };

private: 
  std::shared_ptr<Node> root; 
};

}

template<class mcolor_t>
structure::BinaryTree<mcolor_t>::Node::Node(Node * const par, std::string const & tree, 
    std::unordered_map<std::string, size_t> const & genome_number, 
    std::vector<std::string> const & priority_name) 
: parent(par) 
, left_child (nullptr)
, right_child (nullptr) 
{
  int i = 0; 
  for (i = tree.length() - 1; tree[i] != ':' && tree[i] != ')' && tree[i] != '}' && i >= 0; --i) 
    ; 

  std::string new_tree = tree; 
  if (i > 0) { 
    name = tree.substr(i + 1, tree.length() - i - 1);
    if (tree[i] == ':') { 
      new_tree = tree.substr(0, i);
    } else { 
      new_tree = tree.substr(0, i + 1);
    }  
  }
  

  if (new_tree[0] == '(') {
    //non-trivial tree
    if (new_tree[new_tree.size() - 1] != ')') {
      std::cerr << "ERROR: Bad format input (sub)tree. Check count \')\'" << std::endl;
      exit(3);
    }

    int p = 0;
    for(size_t j = 1; j < new_tree.size() - 1; ++j) {
      if (new_tree[j] == '(' || new_tree[j] == '{') { 
      	++p; 
      } else if (new_tree[j] == ')' || new_tree[j] == '}') {
      	--p;
      } else if (new_tree[j] == ',') { 
      	if (p == 0) {
          left_child = std::unique_ptr<Node>(new Node(this, new_tree.substr(1, j - 1), genome_number, priority_name)); 
      	  right_child = std::unique_ptr<Node>(new Node(this, new_tree.substr(j + 1, new_tree.size() - j - 2), genome_number, priority_name));
      	} 
      } 
      if (p < 0) {
      	std::cerr << "ERROR: Bad format input (sub)tree. Check count \'(\' and \')\'" << std::endl;
      	exit(3);
      }
    }

    if (p != 0) {
      std::cerr << "ERROR: Bad format input (sub)tree. Check count \'(\' and \')\'" << std::endl;
      exit(3);
    }

    if (name.empty()) {
      name = left_child->name + right_child->name;
    } 
    data = mcolor_t(left_child->data, right_child->data, mcolor_t::Union);
  } else {
    if (new_tree[0] == '{' && new_tree[new_tree.size() - 1] == '}') {
      size_t start = 1;
      std::string new_name = "";
      for(size_t j = 1; j < new_tree.size(); ++j) {
        if (new_tree[j] == ',' || new_tree[j] == '}') {
          std::string const & str = new_tree.substr(start, j - start);
          if (genome_number.count(str) == 0) {
          	std::cerr << "ERROR: Unknown genome in (sub)tree: " << str << std::endl;
      	    exit(3);
          }
          data.insert(genome_number.find(str)->second);
          childs.push_back(priority_name[genome_number.find(str)->second]);
          new_name += priority_name[genome_number.find(str)->second];
          start = j + 1;  
        } 
      } 

      if (name.empty()) {
        name = new_name;
      } 
    } else {
      if (genome_number.count(new_tree) == 0) {
        std::cerr << "ERROR: Unknown genome in (sub)tree: " << new_tree << std::endl;
      	exit(3);
      }
      data.insert(genome_number.find(new_tree)->second);
      name = priority_name[genome_number.find(new_tree)->second]; 
    }   
  }
}

template<class mcolor_t>
std::string structure::BinaryTree<mcolor_t>::Node::get_nodes(std::vector<std::string>& info) const  {
  if (!left_child && !right_child) { 
    if (!childs.empty()) { 
      for(auto const & lc : childs) {
        info.push_back("\t\"" + name + "\" -> \"" + lc + "\";");
      }
    } 
    return name; 
  } 

  if (left_child) {
    std::string const & first = left_child->get_nodes(info); 
    info.push_back("\t\"" + name + "\" -> \"" + first + "\";");
  }
  
  if (right_child) {
    std::string const & second = right_child->get_nodes(info); 
    info.push_back("\t\"" + name + "\" -> \"" + second + "\";");
  } 

  return name;
}

template<class mcolor_t>
void structure::BinaryTree<mcolor_t>::Node::get_dicolors(std::set<mcolor_t>& dicolor) const {
  if (left_child) {
    left_child->get_dicolors(dicolor);
  }

  dicolor.insert(data);

  if (right_child) {
    right_child->get_dicolors(dicolor);
  }
} 

template<class mcolor_t>
void structure::BinaryTree<mcolor_t>::Node::get_median_colors(mcolor_t const & complete_color, std::vector<colors_median_t>& medians) const {
  if (left_child) { 
    left_child->get_median_colors(complete_color, medians); 
  }

  if (right_child) { 
    right_child->get_median_colors(complete_color, medians); 
  }

  if (left_child && right_child && this->parent != nullptr) { 
    std::array<mcolor_t, 6> colors;  
    colors[0] = left_child->data;
    colors[1] = right_child->data;
    colors[3] = this->data;
    colors[2] = mcolor_t(complete_color, colors[3], mcolor_t::Difference);
    colors[4] = mcolor_t(colors[2], colors[0], mcolor_t::Union);
    colors[5] = mcolor_t(colors[2], colors[1], mcolor_t::Union);
    medians.push_back(std::set<mcolor_t>(colors.begin(), colors.end()));
  }
  
}
    

template<class mcolor_t>
void structure::BinaryTree<mcolor_t>::Node::get_name_for_colors(std::map<mcolor_t, std::string>& colors) const {
  if (left_child) {
    left_child->get_name_for_colors(colors);
  }

  colors.insert(std::make_pair(data, name));
  
  if (right_child) {
    right_child->get_name_for_colors(colors);
  }
}

template<class mcolor_t>
template<class linearizator_t>
void structure::BinaryTree<mcolor_t>::Node::walk_and_linearizeate(
      linearizator_t const & linearizator, 
      std::map<mcolor_t, typename linearizator_t::partgraph_t> & graphs, 
      std::map<mcolor_t, typename linearizator_t::transform_t> & transformations) const { 
  if (left_child) { 
    left_child->walk_and_linearizeate(linearizator, graphs, transformations); 
  }

  if (right_child) { 
    right_child->walk_and_linearizeate(linearizator, graphs, transformations); 
  }

  if (left_child && right_child && graphs.find(this->data) != graphs.end()) { 
    size_t count_left = linearizator.count_circular_chromosome(graphs[left_child->data]); 
    size_t central = linearizator.count_circular_chromosome(graphs[this->data]); 
    size_t count_right = linearizator.count_circular_chromosome(graphs[right_child->data]); 
    //std::cerr << "Left have " << count_left << std::endl 
    //      << " Central have " << central << std::endl << "Right have " << count_right << std::endl;

    if (count_left == 0 && count_right == 0 && central != 0) { 
      typedef typename linearizator_t::twobreak_t twobreak_t;
      typedef typename linearizator_t::transform_t transform_t;
      typedef typename linearizator_t::partgraph_t partgraph_t;

      //std::cerr << "Start linearizator " << genome_match::mcolor_to_name(this->data) 
      //  << " -> " << genome_match::mcolor_to_name(left_child->data) << std::endl;
      std::pair<transform_t, transform_t> new_history = linearizator.linearizate(graphs[this->data], transformations[left_child->data], graphs[left_child->data]); 

      /*Apply linearization twobreaks*/      
      for (twobreak_t const & twobreak : new_history.first) {
        //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
        //<< " " << twobreak.get_vertex(3) << p" " << std::endl;    
        twobreak.apply_single(graphs[this->data]); 
        //std::cerr << linearizator.count_circular_chromosome(graphs[this->data]) << std::endl;
      }
      
      /*Check that all is good*/
      //std::cerr << "Result genome have " << linearizator.count_circular_chromosome(graphs[this->data]) << std::endl;
      assert(linearizator.count_circular_chromosome(graphs[this->data]) == 0);

      //std::cerr << "Start to check linearizator " << std::endl;
      partgraph_t current = graphs[this->data]; 
      for (twobreak_t const & twobreak : new_history.second) { 
        //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) 
        //<< " " << twobreak.get_vertex(3) << " " << std::endl;    
        
        twobreak.apply_single(current);
      } 
      //std::cerr << "Check that we get good history " << std::endl;
      assert(current == graphs[left_child->data]);

      /*Modify transformation*/
      //std::cerr << "modify transformation" << std::endl;
      for (twobreak_t const & twobreak : new_history.first) { 
        transformations[this->data].push_back(twobreak);
        transformations[right_child->data].push_front(twobreak.inverse());  
      } 
      transformations[left_child->data] = new_history.second;

      //std::cerr << "apply on right child" << std::endl;
      current = graphs[this->data]; 
      for (twobreak_t const & twobreak : transformations[right_child->data]) { 
        twobreak.apply_single(current);
      } 
      assert(current == graphs[right_child->data]);

      if (this->parent->parent != nullptr) {  
        current = graphs[this->parent->data]; 
        for (twobreak_t const & twobreak : transformations[this->data]) { 
          twobreak.apply_single(current);
        } 
        assert(current == graphs[this->data]);
      } 
    }
  } 
}

#endif
