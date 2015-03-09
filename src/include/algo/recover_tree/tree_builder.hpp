#ifndef TREE_BUILDER_HPP__
#define TREE_BUILDER_HPP__

#include <memory>
#include <cassert>
#include "../../structures/tree.hpp"

namespace algo {

  template <class tree_t>
  struct TreeBuilder {
    using node_t = typename tree_t::colored_node_t;
    using mcolor_t = typename node_t::multicolor_t;
    using node_unique_ptr = typename tree_t::node_unique_ptr;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using statistic_t = std::pair<branch_t, size_t>;
    using statistic_vector = std::vector<statistic_t>;

    // Vector should sorted in the descending format
    TreeBuilder(statistic_vector statistics): m_branch_statistics(statistics) {
      validate_statistics();

      if (m_branch_statistics.empty()) {
        ERROR("No branch statistics supplied")
        exit(4);
      }

      auto iter = m_branch_statistics.begin();
      m_root_node = node_unique_ptr(new node_t(mcolor_t(iter->first.first, iter->first.second, mcolor_t::Union)));

      for (;iter != m_branch_statistics.end(); ++iter) {
        //TODO: algo goes here
      }
    }

    node_unique_ptr get_result() {
      return std::move(m_root_node);
    }

  private:
    void validate_statistics() {
      // Maximum size_t value
      size_t previous_value = static_cast<size_t>(-1);
      for (auto& statistic: m_branch_statistics) {
        assert(previous_value >= statistic.second &&
            "Works only on sorted statistics");
        previous_value = statistic.second;
      }
    }

    statistic_vector& m_branch_statistics;
    node_unique_ptr m_root_node;
  };
}

#endif