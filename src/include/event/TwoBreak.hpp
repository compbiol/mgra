#ifndef TWOBREAK_HPP
#define TWOBREAK_HPP

namespace event { 

template<class mcolor_t>
struct TwoBreak {
  typedef typename mcolor_t::citer citer; 
  typedef std::pair<vertex_t, vertex_t> edge_t;
  typedef utility::sym_multihashmap<vertex_t> partgraph_t;
  
  TwoBreak() { 
  } 
  	
  TwoBreak(edge_t const & a1, edge_t const & a2, mcolor_t const & multicolor)
  : m_arcs({a1, a2}) 
  , m_multicolor(multicolor) 
  {
  }

  TwoBreak(vertex_t const & x1, vertex_t const & y1, vertex_t const & x2, vertex_t const & y2, mcolor_t const & multicolor)
  : m_arcs({edge_t(x1, y1), edge_t(x2, y2)})
  , m_multicolor(multicolor) 
  {
  }
 
  inline void change_vertex(size_t index, vertex_t const & v) {
    if (index == 0) {
      m_arcs[0].first = v; 
    } else if (index == 1) {
      m_arcs[0].second = v; 
    } else if (index == 2) {
      m_arcs[1].first = v; 
    } else if (index == 3) {
      m_arcs[1].second = v; 
    } 
  }
  
  inline TwoBreak inverse() const {
    return TwoBreak(m_arcs[0].first, m_arcs[1].first, m_arcs[0].second, m_arcs[1].second, m_multicolor);
  }

  inline edge_t get_arc(size_t index) const { 
    assert(index < 2); 
    return m_arcs[index];
  } 

  inline vertex_t get_vertex(size_t index) const { 
    if (index == 0) {
      return m_arcs[0].first; 
    } else if (index == 1) {
      return m_arcs[0].second; 
    } else if (index == 2) {
      return m_arcs[1].first; 
    } else if (index == 3) {
      return m_arcs[1].second; 
    } 
    assert(false);
  } 

  void apply_single(partgraph_t& local_graph) const;

  TwoBreak get_canonical_twobreak() const;
   
  bool is_independent(TwoBreak const & tested) const;

  bool operator < (TwoBreak const & second) const; 

  bool operator > (TwoBreak const & second) const { 
    return (second < *this);
  }

  DECLARE_GETTER( mcolor_t, m_multicolor, mcolor )
  
  DECLARE_CONST_ITERATOR( citer, m_multicolor, begin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, m_multicolor, end, cend )
  DECLARE_CONST_ITERATOR( citer, m_multicolor, cbegin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, m_multicolor, cend, cend )

private: 
  std::array<edge_t, 2> m_arcs; // (x1,y1) x (x2,y2) = (x1,x2) + (y1,y2)
  mcolor_t m_multicolor; 
};

}

template<class mcolor_t>
void event::TwoBreak<mcolor_t>::apply_single(partgraph_t& local_graph) const {
  for(size_t i = 0; i < 2; ++i) {
    if (m_arcs[i].first != Infty || m_arcs[i].second != Infty) {
      local_graph.erase(m_arcs[i].first, m_arcs[i].second);
    }
  }
  
  if (m_arcs[0].first != Infty || m_arcs[1].first != Infty) {
    local_graph.insert(m_arcs[0].first, m_arcs[1].first);
  }

  if (m_arcs[0].second != Infty || m_arcs[1].second != Infty) {
    local_graph.insert(m_arcs[0].second, m_arcs[1].second);
  }
} 

template<class mcolor_t>
bool event::TwoBreak<mcolor_t>::operator < (TwoBreak const & twobreak) const { 
  auto const less_lambda = [&](vertex_t const & v, vertex_t const & u) -> bool { 
    if (u == Infty) { 
      return true;
    } else if (v == Infty) {
      return false; 
    } else if (std::stoi(v.substr(0, v.length() - 1)) < std::stoi(u.substr(0, u.length() - 1))) { 
      return true;
    } else if (std::stoi(v.substr(0, v.length() - 1)) == std::stoi(u.substr(0, u.length() - 1))) {
      if (*v.crbegin() < *u.crbegin()) {
        return true; 
      } else { 
        return false; 
      }
    } else { 
      return false; 
    }
  };

  auto const equal = [&](vertex_t const & v, vertex_t const & u) -> bool { 
    return !(less_lambda(v, u) && less_lambda(u, v));
  };

  if (less_lambda(m_arcs[0].first, twobreak.m_arcs[0].first)) {
    return true; 
  } else if (equal(m_arcs[0].first, twobreak.m_arcs[0].first)) {
    if (less_lambda(m_arcs[0].second, twobreak.m_arcs[0].second)) {
      return true; 
    } else if (equal(m_arcs[0].second, twobreak.m_arcs[0].second)) {
      if (less_lambda(m_arcs[1].first, twobreak.m_arcs[1].first)) {
        return true; 
      } else if (equal(m_arcs[1].first, twobreak.m_arcs[1].first)) {
        if (less_lambda(m_arcs[1].second, twobreak.m_arcs[1].second)) {
          return true; 
        } else if (equal(m_arcs[1].second, twobreak.m_arcs[1].second)) {
          return (m_multicolor < twobreak.m_multicolor);
        }
      }   
    }
  }

  return false;
}

template<class mcolor_t>
bool event::TwoBreak<mcolor_t>::is_independent(TwoBreak const & tested) const {
  mcolor_t inter_color(m_multicolor, tested.m_multicolor, mcolor_t::Intersection); 

  if (inter_color.empty()) {
    return true; 
  } else {
    auto check_lambda = [&] (size_t ind1, size_t ind2, size_t ind3) -> bool {
      if (tested.m_arcs[ind1] == std::make_pair(get_vertex(ind2), get_vertex(ind3))) { 
        return false;
      } 
      if (tested.m_arcs[ind1] == std::make_pair(get_vertex(ind3), get_vertex(ind2))) { 
        return false;
      }
      return true;
    };

    bool answer = true;
    for (size_t j = 0; j < 2; ++j) { 
      answer = answer && check_lambda(j, 0, 2) && check_lambda(j, 1, 3);
    }
    return answer;
  }

}

template<class mcolor_t>
event::TwoBreak<mcolor_t> event::TwoBreak<mcolor_t>::get_canonical_twobreak() const { 
  edge_t new_arc[2];    
  
  auto const less_lambda = [&](vertex_t const & v, vertex_t const & u) -> bool { 
    if (u == Infty) { 
      return true;
    } else if (v == Infty) {
      return false; 
    } else if (std::stoi(v.substr(0, v.length() - 1)) < std::stoi(u.substr(0, u.length() - 1))) { 
      return true;
    } else if (std::stoi(v.substr(0, v.length() - 1)) == std::stoi(u.substr(0, u.length() - 1))) {
      if (*v.crbegin() < *u.crbegin()) {
        return true; 
      } else { 
        return false; 
      }
    } else { 
      return false; 
    }
  };

  if (!less_lambda(m_arcs[0].first, m_arcs[0].second) && !less_lambda(m_arcs[1].first, m_arcs[1].second)) {
    new_arc[0] = edge_t(m_arcs[0].second, m_arcs[0].first);          
    new_arc[1] = edge_t(m_arcs[1].second, m_arcs[1].first);          
  } else {
    new_arc[0] = m_arcs[0]; 
    new_arc[1] = m_arcs[1];
  }

  if (new_arc[0].first > new_arc[1].first) { 
    std::swap(new_arc[0], new_arc[1]);
  } 

  
  if (!less_lambda(new_arc[0].first, new_arc[1].first)) { 
    std::swap(new_arc[0], new_arc[1]);
  }
  
  return TwoBreak(new_arc[0], new_arc[1], m_multicolor);
}

#endif
