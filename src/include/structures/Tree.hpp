#ifndef TREE_HPP 
#define TREE_HPP

namespace structure { 

template<class mcolor_t>
struct BinaryTree { 

  struct Node { 
    Node() 
    : parent(nullptr) 
    , left_child (nullptr)
    , right_child (nullptr) 
    { 
    }

    Node(Node* const par, std::string const & tree, 
      std::unordered_map<std::string, size_t> const & genome_number, 
      std::vector<std::string> const & priority_name);     

    mcolor_t const & get_data() const { 
      return data;
    }

    std::string const & get_name() const {
      return name;
    }

    std::vector<std::string> const & get_childs() const { 
      return childs;
    }

    Node const * get_parent() const {
      return parent;
    }

    std::unique_ptr<Node> const & get_left_child() const { 
      return left_child;
    }

    std::unique_ptr<Node> const & get_right_child() const { 
      return right_child;
    }

    template<class linearizator_t>
    void walk_and_linearizeate(
      linearizator_t const & linearizator, 
      std::map<mcolor_t, mcolor_t> & parent_colors,
      std::map<mcolor_t, typename linearizator_t::partgraph_t> & graphs, 
      std::map<mcolor_t, typename linearizator_t::transform_t> & transformations) const;

    mcolor_t data;    
  private: 
    std::string name;   
    std::vector<std::string> childs;

    Node* const parent;    
    std::unique_ptr<Node> left_child;    
    std::unique_ptr<Node> right_child; 
  };

  BinaryTree(std::string const & st, std::unordered_map<std::string, size_t> const & genome_number, std::vector<std::string> const & priority_name) 
  : root(new Node(nullptr, st, genome_number, priority_name))
  {
  } 

  std::unique_ptr<Node> const & get_root() const { 
    return root;
  }

  template<class linearizator_t>
  void walk_and_linearizeate(
    linearizator_t const & linearizator, 
    std::map<mcolor_t, mcolor_t> & parent_colors,
    std::map<mcolor_t, typename linearizator_t::partgraph_t> & graphs, 
    std::map<mcolor_t, typename linearizator_t::transform_t> & transformations) const 
  {
    root->walk_and_linearizeate(linearizator, parent_colors, graphs, transformations);

    mcolor_t left = parent_colors[root->get_left_child()->data];
    mcolor_t right = parent_colors[root->get_right_child()->data];

    if (transformations[left].size() < transformations[right].size()) { 
      parent_colors.erase(left); 
      for (auto const & twobreak : transformations[left]) { 
        transformations[right].push_front(twobreak.inverse());
      }
      transformations.erase(left);
    } else { 
      parent_colors.erase(right); 
      for (auto const & twobreak : transformations[right]) { 
        transformations[left].push_front(twobreak.inverse());
      }
      transformations.erase(right);
    }
  }

private: 
  std::unique_ptr<Node> root; 

private: 
  DECL_LOGGER("BinaryTree") 
};

template<class mcolor_t>
BinaryTree<mcolor_t>::Node::Node(Node * const par, std::string const & tree, 
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
      ERROR("Bad format input (sub)tree. Check count \')\'")
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
      	ERROR("Bad format input (sub)tree. Check count \'(\' and \')\'")
      	exit(3);
      }
    }

    if (p != 0) {
      ERROR("Bad format input (sub)tree. Check count \'(\' and \')\'")
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
          	ERROR("Unknown genome in (sub)tree: " << str)
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
        ERROR("Unknown genome in (sub)tree: " << new_tree)
      	exit(3);
      }
      data.insert(genome_number.find(new_tree)->second);
      name = priority_name[genome_number.find(new_tree)->second]; 
    }   
  }
}

}

/*ERROR, FIXME, think about edge from root in this function*/
template<class mcolor_t>
template<class linearizator_t>
void structure::BinaryTree<mcolor_t>::Node::walk_and_linearizeate(
      linearizator_t const & linearizator, 
      std::map<mcolor_t, mcolor_t> & parent_colors,
      std::map<mcolor_t, typename linearizator_t::partgraph_t> & graphs, 
      std::map<mcolor_t, typename linearizator_t::transform_t> & transformations) const { 
  if (this->parent != nullptr) {  
    if (this->parent->parent != nullptr) { 
      parent_colors.insert(std::make_pair(this->data, this->parent->data));
    } else { 
      if (this->parent->left_child.get() == this) { 
        parent_colors.insert(std::make_pair(this->data, this->parent->right_child->data));
      } else if (this->parent->right_child.get() == this) { 
        parent_colors.insert(std::make_pair(this->data, this->parent->left_child->data));
      }
    }
  } 

  if (left_child) { 
    left_child->walk_and_linearizeate(linearizator, parent_colors, graphs, transformations); 
  }

  if (right_child) { 
    right_child->walk_and_linearizeate(linearizator, parent_colors, graphs, transformations); 
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
