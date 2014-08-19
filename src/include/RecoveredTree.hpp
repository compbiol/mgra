#ifndef RECOVERED_TREE_HPP
#define RECOVERED_TREE_HPP

template<mcolor_t>
struct RecoveredTree { 
  typedef std::tuple<mcolor_t, mcolor_t, size_t> branch_t;
  typedef typename structure::BinaryTree<mcolor_t> tree_t;

  RecoveredTree(std::vector<branch_t> const & colors_information, std::vector<tree_t> const & input_subtrees) 
  : m_colors_information(colors_information)
  , m_input_subtrees(input_subtrees)
  {
  }

  
private:
  std::vector<branch_t> m_colors_information;
  std::list<tree_t> m_input_subtrees;
  std::vector<std::pair<tree_t, size_t> > m_output_tree; 
};

#endif 