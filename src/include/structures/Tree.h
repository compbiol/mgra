#ifndef TREE_H_ 
#define TREE_H_

namespace structure { 

template<class mcolor_t>
struct BinaryTree { 
  struct Node { 
    std::string name;   
    std::vector<std::string> childs;
    mcolor_t data;
    std::shared_ptr<Node> left_child;    
    std::shared_ptr<Node> right_child; 

    Node () 
    : left_child (nullptr)
    , right_child (nullptr) 
    { 
    }

    Node(std::string const & tree, std::unordered_map<std::string, size_t> const & genome_number, std::vector<std::string> const & priority_name);     

    std::string get_nodes(std::vector<std::string>& info) const;

    void get_dicolors(std::set<mcolor_t>& dicolor) const; 
    void get_name_for_colors(std::map<mcolor_t, std::string>& colors) const;
  };

  BinaryTree(std::string const & st, std::unordered_map<std::string, size_t> const & genome_number, std::vector<std::string> const & priority_name) 
  : root(st, genome_number, priority_name)
  { 
  } 

  void get_nodes(std::vector<std::string>& info) const {
    root.get_nodes(info);	
  }

  std::set<mcolor_t> build_vec_T_consistent_colors() const {
    std::set<mcolor_t> dicolor;
    root.get_dicolors(dicolor);
    return dicolor;	
  }

  std::map<mcolor_t, std::string> get_name_for_colors() const {
    std::map<mcolor_t, std::string> colors;
    root.get_name_for_colors(colors);
    return colors;	
  }

private: 
  Node root; 
};

}

template<class mcolor_t>
structure::BinaryTree<mcolor_t>::Node::Node(std::string const & tree, std::unordered_map<std::string, size_t> const & genome_number, std::vector<std::string> const & priority_name) {
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
          left_child = std::make_shared<Node>(Node(new_tree.substr(1, j - 1), genome_number, priority_name)); 
      	  right_child = std::make_shared<Node>(Node(new_tree.substr(j + 1, new_tree.size() - j - 2), genome_number, priority_name));
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
    left_child = nullptr;
    right_child = nullptr;
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
void structure::BinaryTree<mcolor_t>::Node::get_name_for_colors(std::map<mcolor_t, std::string>& colors) const {
  if (left_child) {
    left_child->get_name_for_colors(colors);
  }

  colors.insert(std::make_pair(data, name));
  
  if (right_child) {
    right_child->get_name_for_colors(colors);
  }
}
 
#endif
