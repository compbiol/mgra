#ifndef INSDEL_H_
#define INSDEL_H_

namespace event { 

template<class mcolor_t>
struct InsDel {
  typedef std::pair<vertex_t, vertex_t> arc_t;
  typedef typename mcolor_t::citer citer; 
	
  InsDel(arc_t const & edge, mcolor_t const & multicolor)
  : m_edge(edge) 
  , m_multicolor(multicolor) 
  {
  }

  InsDel(vertex_t const & x1, vertex_t const & y1, mcolor_t const & multicolor)
  : m_edge(arc_t(x1, y1)) 
  , m_multicolor(multicolor) 
  {
  }
    
  inline arc_t get_edge() const { 
    return m_edge;
  } 

  inline mcolor_t get_mcolor() const { 
    return m_multicolor; 
  } 

  DECLARE_CONST_ITERATOR( citer, m_multicolor, begin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, m_multicolor, end, cend )
  DECLARE_CONST_ITERATOR( citer, m_multicolor, cbegin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, m_multicolor, cend, cend )

private: 
  arc_t m_edge; // (x1, y1)
  mcolor_t m_multicolor; 
};

}
#endif

