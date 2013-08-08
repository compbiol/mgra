#ifndef TREE_H_ 
#define TREE_H_

#include "mcolor.h"

template<class type_data>
struct BinaryTree { 
  struct Node { 
    //std::vector<std::string> data;     
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

    Mcolor get_dicolors(std::set<Mcolor>& dicolor) const; 
  };

  BinaryTree(const std::string& st, const std::unordered_map<std::string, size_t>& genome_number) 
  : root(st, genome_number)
  { 
  } 

  template<class cfg_t>
  void get_nodes(std::vector<std::string>& info, const cfg_t& cfg) const {
	root.get_nodes(info, cfg);	
  }

  void build_vec_T_consistent_colors(std::set<Mcolor>& dicolor) const {
	Mcolor color = root.get_dicolors(dicolor);	
	dicolor.insert(color); 
  }

private: 
  Node root; 
};

template<class type_data>
BinaryTree<type_data>::Node::Node(const std::string& tree, const std::unordered_map<std::string, size_t>& genome_number) {
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
std::string/*std::set<std::string>*/ BinaryTree<type_data>::Node::get_nodes(std::vector<std::string>& info, const cfg_t& cfg) const  {
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
Mcolor BinaryTree<type_data>::Node::get_dicolors(std::set<Mcolor>& dicolor) const {
	if (!left_child && !right_child) { 
		Mcolor temp; 
		for(auto it = data.cbegin(); it != data.cend(); ++it) {	
			temp.insert(*it); 
    		}
		return temp;
	} 

	Mcolor first;  
	if (left_child) {
		first = left_child->get_dicolors(dicolor);
		dicolor.insert(first);
	}

	Mcolor second; 
	if (right_child) {
		second = right_child->get_dicolors(dicolor);
		dicolor.insert(second);
	}

	Mcolor result(first, second, Mcolor::Union);

	return result;
} 
#endif
