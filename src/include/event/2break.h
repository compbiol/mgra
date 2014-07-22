#ifndef TWOBREAK_H_
#define TWOBREAK_H_

#include "Event.h"


namespace event { 

template<class mcolor_t>
struct TwoBreak : public event::Event {
  typedef std::pair<vertex_t, vertex_t> arc_t;
  typedef typename mcolor_t::citer citer; 

  TwoBreak() { 
  } 
  	
  TwoBreak(arc_t const & a1, arc_t const & a2, mcolor_t const & multicolor)
  : m_arcs({a1, a2}) 
  , m_multicolor(multicolor) 
  {
  }

  TwoBreak(vertex_t const & x1, vertex_t const & y1, vertex_t const & x2, vertex_t const & y2, mcolor_t const & multicolor)
  : m_arcs({arc_t(x1, y1), arc_t(x2, y2)})
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

  inline arc_t get_arc(size_t index) const { 
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

  inline mcolor_t get_mcolor() const { 
    return m_multicolor; 
  } 

  void apply_single(partgraph_t& local_graph) const;

  DECLARE_CONST_ITERATOR( citer, m_multicolor, begin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, m_multicolor, end, cend )
  DECLARE_CONST_ITERATOR( citer, m_multicolor, cbegin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, m_multicolor, cend, cend )

private: 
  //std::array<vertex_t, 4> arcs;
  arc_t m_arcs[2]; // (x1,y1) x (x2,y2) = (x1,x2) + (y1,y2)
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

#endif

