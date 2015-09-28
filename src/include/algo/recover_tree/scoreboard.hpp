//
// Created by Nikita Kartashov on 19/03/2015.
//

#ifndef MGRA_SCOREBOARD_HPP_
#define MGRA_SCOREBOARD_HPP_

namespace algo {
  template <class value_t>
  struct Scoreboard {
    using scored_position_t = std::pair<size_t, value_t>;

    Scoreboard(size_t max_scores) {
      assert(max_scores != 0);
      m_scores.resize(max_scores);
    }

    scored_position_t const& top() const {
      return m_scores[0];
    }

    scored_position_t const& operator[](size_t index) const {
      return m_scores[index];
    }

    size_t min_score() const {
      return m_scores.back().first;
    }

    size_t max_score() const {
      return m_scores.front().first;
    }

    bool update(scored_position_t scored_position) {
      if (scored_position.first <= min_score()) {
        return false;
      }
      if (scored_position.first > max_score()) {
        std::rotate(std::begin(m_scores), std::begin(m_scores) + 1, std::end(m_scores));
        m_scores[0] = scored_position;
        return true;
      }

      // Lower bound cannot point to the last element, it would have been rejected by 1st if clause
      auto lower_bound = std::lower_bound(m_scores.rbegin(), m_scores.rend(), scored_position).base();
      std::rotate(lower_bound, lower_bound + 1, std::end(m_scores));

      return true;
    }

  private:
    std::vector<scored_position_t> m_scores;
  };
}

#endif //_MGRA_SCOREBOARD_HPP_
