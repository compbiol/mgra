//
// Created by Nikita Kartashov on 18/03/2015.
//

#ifndef MGRA_MCOLOR_INFO_HPP_
#define MGRA_MCOLOR_INFO_HPP_

namespace structure {
  template <class mcolor_t>
  struct McolorInfo {
    /**
    * Position in color column
    * first: level index (always > 0)
    * second: position on level
    */
    using position_t = std::pair<size_t, size_t>;
    /**
    * Two positions subtree colors, making a tree
    */
    using subtree_positions_t = std::pair<position_t, position_t>;

    McolorInfo() {
    }

    McolorInfo(size_t separate_score) :
        m_separate_score(separate_score),
        m_whole_tree_score(m_separate_score),
        m_subtree_positions(subtree_positions_t()) {
      m_separate_score = separate_score;
    }

    McolorInfo(
        size_t separate_score,
        size_t whole_tree_score,
        subtree_positions_t positions) :
        m_separate_score(separate_score),
        m_whole_tree_score(whole_tree_score),
        m_subtree_positions(positions) {
    }

    /**
    * Updates info with the score & positions of a better tree, otherwise nothing
    * @return true if updated
    */
    bool update_if_better(McolorInfo const& info) {
      if (m_whole_tree_score >= info.m_whole_tree_score) {
        return false;
      }
      m_whole_tree_score = info.m_whole_tree_score;
      m_subtree_positions = info.m_subtree_positions;
      return true;
    }

    McolorInfo make_from_subtree_info(
        McolorInfo const& left,
        McolorInfo const& right,
        position_t const& left_position,
        position_t const& right_position) const {
      auto new_score = m_separate_score + left.m_whole_tree_score + right.m_whole_tree_score;
      return McolorInfo(m_separate_score, new_score, std::make_pair(left_position, right_position));
    }

    size_t get_separate_score() const {
      return m_separate_score;
    }

    size_t get_whole_tree_score() const {
      return m_whole_tree_score;
    }

    subtree_positions_t const& get_subtree_positions() const {
      return m_subtree_positions;
    }

    bool has_subtrees() const {
      return m_separate_score < m_whole_tree_score;
    }

  private:
    size_t m_separate_score;
    size_t m_whole_tree_score;
    subtree_positions_t m_subtree_positions;
  };
}

#endif //_MGRA_MCOLOR_INFO_HPP_
