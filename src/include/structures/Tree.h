#ifndef TREE_H_ 
#define TREE_H_

#include "mcolor.h"

namespace structure { 

template<class type_data>
struct BinaryTree { 
  typedef structure::Mcolor mcolor_t;

  struct Node { 
    std::string name;   
    std::vector<size_t> data;
    std::shared_ptr<Node> left_child;    
    std::shared_ptr<Node> right_child; 

    Node(const std::string& tree, const std::unordered_map<std::string, size_t>& genome_number);     

    Node () 
    : left_child (nullptr)
    , right_child (nullptr) 
    { 
    }

    template<class cofg_t>
    std::string/*std::set<std::string>*/ get_nodes(std::vector<std::string>& info, const cofg_t& cfg) const;

    mcolor_t get_dicolors(std::set<mcolor_t>& dicolor) const; 
  };

  BinaryTree(const std::string& st, const std::unordered_map<std::string, size_t>& genome_number) 
  : root(st, genome_number)
  { 
  } 

  template<class cfg_t>
  void get_nodes(std::vector<std::string>& info, const cfg_t& cfg) const {
	root.get_nodes(info, cfg);	
  }

  void build_vec_T_consistent_colors(std::set<mcolor_t>& dicolor) const {
	const mcolor_t& color = root.get_dicolors(dicolor);	
	dicolor.insert(color); 
  }

private: 
  Node root; 
};

}

template<class type_data>
structure::BinaryTree<type_data>::Node::Node(const std::string& tree, const std::unordered_map<std::string, size_t>& genome_number) {
  int i = 0; 
  for (i = tree.length() - 1; tree[i] != ':' && tree[i] != ')' && i >= 0; --i) 
    ; 

  std::string new_tree = tree; 
  if (i > 0 && tree[i] == ':') {
    name = tree.substr(i + 1, tree.length() - i - 1);
    new_tree = tree.substr(0, i); 
  } 

  if (tree[0] == '(') {
    //non-trivial tree
    if (tree[tree.size() - 1] != ')') {
      std::cerr << "ERROR: Malformed input (sub)tree 1" << std::endl;
      exit(3);
    }

    int p = 0;
    for(size_t j = 1; j < tree.size() - 1; ++j) {
      if (tree[j] == '(') { 
	++p; 
      } else if (tree[j] == ')') {
	--p;
      } else if (tree[j] == ',') { 
	if (p == 0) {
	  left_child = std::make_shared<Node>(Node(tree.substr(1, j - 1), genome_number)); 
	  right_child = std::make_shared<Node>(Node(tree.substr(j + 1, tree.size() - j - 2), genome_number));
	} 
      } 
      if (p < 0) {
	std::cerr << "ERROR: Malformed input (sub)tree 2" << std::endl;
	exit(3);
      }
    }
    if (p != 0) {
      std::cerr << "ERROR: Malformed input (sub)tree 3" << std::endl;
      exit(3);
    }
  } else {
    /*
    if (new_tree[0] == '{' && new_tree[new_tree.size() - 1] == '}') {
      size_t start = 1;
      for(size_t j = 1; j < new_tree.size() - 1; ++j) {
        if (new_tree[j] == ',') {
          std::string str = new_tree.substr(start, new_tree.size() - start - 1);
          data.push_back(genome_number.find(str)->second);
          start = j + 1;  
        } 
      } 
    } else {
      if (genome_number.count(new_tree) == 0) {
	data.push_back(genome_number.find(new_tree)->second);
      } 
    } 
    */
      //single node
      for(size_t j = 0; j < tree.size(); ++j) {
        std::string c = tree.substr(j, 1);
        if (genome_number.count(c) == 0) {
	  std::cerr << "ERROR: Unknown genome in (sub)tree: " << tree << std::endl;
	  exit(3);
        }
        data.push_back(genome_number.find(c)->second);
      }
   
    left_child = nullptr;
    right_child = nullptr;
  }
}

template<class type_data>
template<class cfg_t>
std::string/*std::set<std::string>*/ structure::BinaryTree<type_data>::Node::get_nodes(std::vector<std::string>& info, const cfg_t& cfg) const  {
	if (!left_child && !right_child) { 
		//std::set<std::string> lists;
		std::string temp = ""; 
		for(auto it = data.cbegin(); it != data.cend(); ++it) {
			temp += (cfg.get_priority_name(*it));//lists.insert(*it);
    		}
		return temp;//return lists;
	} 

	std::string first = "";  
	//std::set<std::string> list_first; 
	if (left_child) {
		first = left_child->get_nodes(info, cfg);
		//list_first = left_child->get_nodes(info);
	
		/*for (auto it = list_first.cbegin(); it != list_first.cend(); ++it) {
			first += (*it);
		}*/
	}

	std::string second = ""; 
	//std::set<std::string> list_second; 
	if (right_child) {
		second = right_child->get_nodes(info, cfg);
		//list_second = right_child->get_nodes(info);
		/*for (auto it = list_second.cbegin(); it != list_second.cend(); ++it) {
			second += (*it);
		}*/
	}

	/*std::set<std::string> list_result;
	std::set_union(list_first.cbegin(), list_first.cend(), list_second.cbegin(), list_second.cend(), std::inserter(list_result, list_result.begin()));*/

	std::string result = first + second; ///"";
	//for (auto it = list_result.cbegin(); it != list_result.cend(); ++it) {
	//	result += (*it);
	//}

	info.push_back("\t\"" + result + "\"\t->\t\"" + first + "\";");
	info.push_back("\t\"" + result + "\"\t->\t\"" + second + "\";");

	//return list_result;
	return result;
}

template<class type_data>
structure::Mcolor structure::BinaryTree<type_data>::Node::get_dicolors(std::set<mcolor_t>& dicolor) const {
	if (!left_child && !right_child) { 
		mcolor_t temp; 
		for(auto it = data.cbegin(); it != data.cend(); ++it) {	
			temp.insert(*it); 
    		}
		return temp;
	} 

	mcolor_t first;  
	if (left_child) {
		first = left_child->get_dicolors(dicolor);
		dicolor.insert(first);
	}

	mcolor_t second; 
	if (right_child) {
		second = right_child->get_dicolors(dicolor);
		dicolor.insert(second);
	}

	mcolor_t result(first, second, mcolor_t::Union);

	return result;
} 

#endif
