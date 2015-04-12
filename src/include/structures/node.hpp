#ifndef NODE_HPP_
#define NODE_HPP_

#include <iostream>
#include <memory>
#include <vector>
#include <unordered_map>

#include "logger/logger.hpp"

namespace structure {

  template <class mcolor_t>
  struct Node {
    using multicolor_t = mcolor_t;
    using node_t = Node<mcolor_t>;
    using node_ptr = std::shared_ptr<node_t>;
    using node_ptr_ref = node_ptr&;
    using node_unique_ptr = std::unique_ptr<node_t>;
    using node_const_ptr = node_ptr const&;

    Node(mcolor_t color) : data(color),
                           m_complete(color.size() == 1),
                           parent(nullptr),
                           left_child(nullptr),
                           right_child(nullptr) {
    }

    Node(Node* const par, std::string const& tree,
         std::unordered_map<std::string, size_t> const& genome_number,
         std::vector<std::string> const& priority_name);

    mcolor_t const& get_data() const {
      return data;
    }

    std::string const& get_name() const {
      return name;
    }

    std::vector<std::string> const& get_children() const {
      return children;
    }

    Node const* get_parent() const {
      return parent;
    }

    node_const_ptr get_left_child() const {
      return left_child;
    }

    node_const_ptr get_right_child() const {
      return right_child;
    }

    node_ptr_ref get_left_child() {
      return left_child;
    }

    node_ptr_ref get_right_child() {
      return right_child;
    }

    void set_parent(Node* new_parent) {
      parent = new_parent;
    }

    void set_name(std::string new_name) {
      name = new_name;
    }


    void set_left_child(node_ptr node) {
      left_child = std::move(node);
    }

    void set_right_child(node_ptr node) {
      right_child = std::move(node);
    }

    bool has_left_child() const {
      return static_cast<bool>(left_child);
    }

    bool has_right_child() const {
      return static_cast<bool>(right_child);
    }

    size_t get_size() const {
      return data.size();
    }

    bool is_simple() const {
      return get_size() == 1;
    }

    bool is_leaf() const {
      return !has_left_child() && !has_right_child();
    }

    bool is_root() const {
      return parent == nullptr;
    }

    bool is_complete() {
      if (m_complete) {
        return m_complete;
      };
      if (has_left_child()) {
        m_complete = left_child->is_complete();
      }
      if (!m_complete) {
        return m_complete;
      }
      if (has_right_child()) {
        m_complete = m_complete && right_child->is_complete();
      }
      return m_complete;
    }

    void set_complete() {
      m_complete = true;
    }

  private:
    mcolor_t data;
    std::string name;
    bool is_whole_duplication;
    std::vector<std::string> children;

    /**
    * If the tree downwards has all the leaf nodes with only one color
    */
    bool m_complete;
    Node* parent;
    node_ptr left_child;
    node_ptr right_child;
  };

  template <class mcolor_t>
  Node<mcolor_t>::Node(Node* const par, std::string const& tree,
                       std::unordered_map<std::string, size_t> const& genome_number,
                       std::vector<std::string> const& priority_name)
      : parent(par), left_child(nullptr), right_child(nullptr) {
    size_t i = 0;
    for (i = tree.length() - 1; tree[i] != ':' && tree[i] != ')' && tree[i] != '}' && i >= 0; --i);

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
      for (size_t j = 1; j < new_tree.size() - 1; ++j) {
        if (new_tree[j] == '(' || new_tree[j] == '{') {
          ++p;
        } else if (new_tree[j] == ')' || new_tree[j] == '}') {
          --p;
        } else if (new_tree[j] == ',') {
          if (p == 0) {
            left_child = std::make_shared<node_t>(this, new_tree.substr(1, j - 1), genome_number, priority_name);
            right_child = std::make_shared<node_t>(this, new_tree.substr(j + 1, new_tree.size() - j - 2), genome_number,
                                                   priority_name);
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
        for (size_t j = 1; j < new_tree.size(); ++j) {
          if (new_tree[j] == ',' || new_tree[j] == '}') {
            std::string const& str = new_tree.substr(start, j - start);
            if (genome_number.count(str) == 0) {
              ERROR("Unknown genome in (sub)tree: " << str)
              exit(3);
            }
            data.insert(genome_number.find(str)->second);
            children.push_back(priority_name[genome_number.find(str)->second]);
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

#endif