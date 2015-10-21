#ifndef MCOLOR_HPP
#define MCOLOR_HPP

#include "utility/hash.hpp"

namespace structure {

struct Mcolor {

  enum Construct {
    Difference, Union, Intersection
  };

  using map_t = std::map<size_t, size_t>;
  using color_pair_t = std::pair<size_t, size_t>;
  using citer = map_t::const_iterator;
  using iter = map_t::iterator;
  using branch_t = std::pair<Mcolor, Mcolor>;
  using color_vector_t = std::vector<Mcolor>;
  
  Mcolor() = default; 
  
  explicit Mcolor(size_t color, size_t multiplicity = 1) {
    main_color.insert(std::make_pair(color, multiplicity));
  }

  Mcolor(Mcolor const& first, Mcolor const& second, Construct const& what) {
    switch (what) {
      case Difference:
        set_difference(first, second);
        break;
      case Union:
        set_union(first, second);
        break;
      case Intersection:
        set_intersection(first, second);
        break;
    }
  }

  static Mcolor color_union(color_vector_t const& colors) {
    Mcolor result_color;
    for (auto& color: colors) {
      result_color.add_union(color);
    }
    return result_color;
  }

  color_vector_t break_into_parts() const {
    color_vector_t result;
    std::transform(std::begin(main_color), std::end(main_color), std::back_inserter(result),
                   [](color_pair_t pair) {
                     return Mcolor(pair.first, pair.second);
                   });
    return result;
  }


  void insert(size_t i) {
    if (main_color.find(i) == main_color.end()) {
      main_color.insert(std::make_pair(i, 1));
    } else {
      main_color[i] += 1;
    }
  }

  void clear() {
    main_color.clear();
  }

  bool is_one_to_one_match() const {
    for (auto const& col : main_color) {
      if (col.second != 1) {
        return false;
      }
    }
    return true;
  }

  size_t how_much_includes(Mcolor const& second) const {
    size_t answer = 0;
    Mcolor current = *this;
    while (current.includes(second)) {
      ++answer;
      current = Mcolor(current, second, Mcolor::Difference);
    }
    return answer;
  }

  /**
  * Tells if this color is a superset of a second one
  */
  bool includes(Mcolor const& second) const {
    auto first1 = main_color.cbegin();
    auto first2 = second.cbegin();

    for (; first2 != second.cend(); ++first1) {
      if (first1 == main_color.cend()) {
        return false;
      }

      if (first2->first < first1->first) {
        return false;
      }

      if (first1->first == first2->first) {
        if (first1->second < first2->second) {
          return false;
        } else {
          ++first2;
        }
      }
    }
    return true;
  }

  bool defined(size_t i) const {
    return (main_color.count(i) != 0);
  }

  bool empty() const {
    return main_color.empty();
  }

  size_t size() const {
    size_t size = 0;
    for (auto const& col : main_color) {
      size += col.second;
    }
    return size;
  }

  static branch_t pack(branch_t const& colors) {
    return colors.first.pack(colors.second);
  }

  branch_t pack(Mcolor const& right) const {
    return std::make_pair(std::max(*this, right), std::min(*this, right));
  };

  /**
  * Gets a compliment color, given a superset and then packs it so, the first pair color is larger
  */
  branch_t packed_compliment(Mcolor const& supercolor) const {
    assert(supercolor.includes(*this));
    Mcolor complement(supercolor, *this, Difference);
    return pack(complement);
  }

  size_t make_hash() const {
    size_t seed = 0;

    for (auto const& value: main_color) {
      util::hash_combine(seed, value.first);
      util::hash_combine(seed, value.second);
    }
    return seed;
  }

  bool operator>(Mcolor const& C) const {
    return (main_color > C.main_color);
  }

  bool operator<(Mcolor const& C) const {
    return (main_color < C.main_color);
  }

  bool operator==(Mcolor const& C) const {
    return (main_color == C.main_color);
  }

  bool operator!=(Mcolor const& C) const {
    return (main_color != C.main_color);
  }

  DECLARE_CONST_ITERATOR(citer, main_color, begin, cbegin)

  DECLARE_CONST_ITERATOR(citer, main_color, end, cend)

  DECLARE_CONST_ITERATOR(citer, main_color, cbegin, cbegin)

  DECLARE_CONST_ITERATOR(citer, main_color, cend, cend)

private:
  void set_difference(Mcolor const& first, Mcolor const& second) {
    auto first1 = first.cbegin();
    auto first2 = second.cbegin();
    auto result = std::inserter(main_color, main_color.begin());

    while (first1 != first.cend() && first2 != second.cend()) {
      if (first1->first == first2->first) {
        if (first1->second > first2->second) {
          *result = std::make_pair(first1->first, first1->second - first2->second);
        }
        ++first1;
        ++first2;
      } else if (*first1 < *first2) {
        *result = *first1;
        ++result;
        ++first1;
      } else if (*first2 < *first1) {
        ++first2;
      }
    }

    for (; first1 != first.cend(); ++first1, ++result) {
      *result = *first1;
    }
  }

  void add_union(Mcolor const& color) {
    auto this_iter = color.cbegin();
    auto color_iter = color.cbegin();
    while (true) {
      if (color_iter == color.cend()) {
        break;
      }

      if (this_iter == main_color.cend()) {
        while (color_iter != color.cend()) {
          main_color.insert(*this_iter++);
        }
      }

      if (this_iter->first == color_iter->first) {
        main_color[this_iter->first] += color_iter->second;
        ++this_iter;
        ++color_iter;
      } else if (this_iter->first > color_iter->first) {
        main_color.insert(*color_iter++);
      } else if (this_iter->first < color_iter->first) {
        ++this_iter;
      }
    }
  }

  void set_union(Mcolor const& first, Mcolor const& second) {
    auto first1 = first.cbegin();
    auto first2 = second.cbegin();
    auto result = std::inserter(main_color, main_color.begin());

    while (true) {
      if (first1 == first.cend()) {
        for (; first2 != second.cend(); ++first2, ++result) {
          *result = *first2;
        }
        break;
      }

      if (first2 == second.cend()) {
        for (; first1 != first.cend(); ++first1, ++result) {
          *result = *first1;
        }
        break;
      }

      if (first1->first == first2->first) {
        *result = std::make_pair(first1->first, first1->second + first2->second);
        ++first1;
        ++first2;
      } else if (*first1 < *first2) {
        *result = *first1;
        ++first1;
      } else if (*first2 < *first1) {
        *result = *first2;
        ++first2;
      }
      ++result;
    }
  }

  void set_intersection(Mcolor const& first, Mcolor const& second) {
    auto first1 = first.cbegin();
    auto first2 = second.cbegin();
    auto result = std::inserter(main_color, main_color.begin());

    while (first1 != first.cend() && first2 != second.cend()) {
      if (first1->first == first2->first) {
        *result = std::make_pair(first1->first, std::min(first1->second, first2->second));
        ++result;
        ++first1;
        ++first2;
      } else if (*first1 < *first2) {
        ++first1;
      } else if (*first2 < *first1) {
        ++first2;
      }
    }
  }

  map_t main_color;
};

}

#endif  /* end of include guard: MCOLOR_HPP */
