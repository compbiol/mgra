#ifndef TREE_H_ 
#define TREE_H_

#include <vector>
#include <string>
#include <memory>

#include "mcolor.h"
#include "pconf.h"

template<class type_data>
struct BinaryTree { 
private: 
  struct Node { 
    std::vector<type_data> data; 
    std::shared_ptr<Node> left_child;    
    std::shared_ptr<Node> right_child; 
    
    std::string to_string() const; 
    Mcolor to_color(std::vector<Mcolor>& dicolors) const;
  }; 

  Node root; 

public: 
  BinaryTree(const std::string& st) { 
    std::shared_ptr<Node> my = construct_tree(st);
    root = *my;
  } 

  std::shared_ptr<Node> construct_tree(const std::string& st); 

  std::vector<Mcolor> get_dicolors() { 
    std::vector<Mcolor> ans; 
    root.to_color(ans);
    return ans;
  } 

  std::string print() { 
    return root.to_string();
  } 

};

template<class type_data>
std::shared_ptr<typename BinaryTree<type_data>::Node> BinaryTree<type_data>::construct_tree(const std::string& st) { 
  std::shared_ptr<Node> current(new Node);

  if (st[0] == '(' && st[st.size() - 1] == ')') {  
    std::string temp = ""; 
    size_t count_brackets = 0; 

    for(size_t i = 1; i < st.size() - 1; ++i) { 
      if (st[i] == '(' || st[i] == '{') { 
	++count_brackets;
      } else if (st[i] == ')' || st[i] == '}') {  
	--count_brackets;
      } else if (st[i] == ',') { 
	if (count_brackets == 0) { 
	  current->left_child = construct_tree(st.substr(1, i - 1));
	  current->right_child = construct_tree(st.substr(i + 1, st.length() - i - 2));
	} 
      } 

      if (count_brackets < 0) { 
	std::cerr << "ERROR: Malformed input (sub)tree 3" << std::endl;
	exit(1);	
      } 
    } 
  } else if (st[0] == '{' && st[st.size() - 1] == '}') { 
    std::string temp = ""; 
    for(size_t i = 1; i < st.size() - 1; ++i) { 
      if (st[i] == ',') {
	current->data.push_back(temp);
	temp.clear();
      } else if (check_symbol(st[i])) { 
	temp += st[i]; 
      } else {  
	std::cerr << "ERROR: Malformed input (sub)tree 4" << std::endl;
	exit(1);
      } 
    } 
    current->data.push_back(temp);
  } else if (!st.empty()) { 
      current->data.push_back(st);
  } 

  return current; 
} 

template<class type_data>
Mcolor BinaryTree<type_data>::Node::to_color(std::vector<Mcolor>& dicolors) const {
  if (!left_child && ! right_child) { 	
    Mcolor ans; 
    for (auto it = data.vbegin(); it != data.cend(); ++it) { 
      if (ProblemInstance::member_name(data)) { 
	ans.insert(*it); 
      } else { 
	std::cerr << "ERROR: Unknown genome in (sub)tree: " << *it << std::endl; 
	exit(1);
      } 
    } 
    dicolors.push_back(ans); 
    return ans;
  } 

  Mcolor left;
  Mcolor right;	
  if (left_child) { 
    left = left_child->to_color(dicolors);
  }
  
  if (right_child) { 
    right = right_child->to_color(dicolors);
  }      

  Mcolor ans(left, right, Mcolor::Union);
  dicolors.push_back(ans); 
  return ans; 
} 

template<class type_data>
std::string BinaryTree<type_data>::Node::to_string() const { 
  std::string answer = "("; 

  if (left_child) {
    answer += left_child->to_string();
  } 

  if (data.size() == 1) { 
    answer += (data[0] + ',');
  } else if (!data.empty()) { 
    answer += '{';
    for (auto it = data.cbegin(); it != data.cend(); ++it) { 
      answer += (*it + ','); 
    } 
    answer += '}'; 
  }

  if (right_child) {
    answer += right_child->to_string();
  } 

  return answer + ')'; 
}  

#endif
