#ifndef MULARCS_H_
#define MULARCS_H_

namespace structure { 

template<class mcolor_t>
struct Mularcs { 
  typedef std::multimap<vertex_t, mcolor_t> map_t; 
  typedef typename map_t::const_iterator citer; 
  typedef typename map_t::const_reverse_iterator criter;
  typedef typename map_t::iterator iter;    

  inline void insert(vertex_t const & v, mcolor_t const & mc) { 
    m_mularcs.insert(std::make_pair(v, mc));	
  } 

  inline void insert(vertex_t const & v, size_t i) {
    if (m_mularcs.find(v) != m_mularcs.end()) { 
      m_mularcs.find(v)->second.insert(i);
    } else { 
      m_mularcs.insert(std::make_pair(v, mcolor_t(i)));
    }
  }
 
  inline void erase(vertex_t const & v) { 
    m_mularcs.erase(v);
  }
 
  inline void erase(vertex_t const & v, mcolor_t const & color) { 
    for (auto iter = m_mularcs.begin(); iter != m_mularcs.end();) { 
      if (iter->first == v && iter->second == color) { 
        m_mularcs.erase(iter++);    
      } else { 
        ++iter;
      }
    }    
  }

  inline vertex_t get_vertex(mcolor_t const & color) const { //FIXME: MAYBE OPTIONAL
    vertex_t v = "";
    for (auto const & arc : m_mularcs) {
      if (arc.second == color) { 
        v = arc.first;
        break;
      }
    }
    return v;   
  } 
  
  inline std::pair<citer, citer> equal_range(vertex_t const & v) const {
    return m_mularcs.equal_range(v);
  } 
 
  inline mcolor_t union_multicolors() const {
    mcolor_t un; 
    for (auto const &arc : m_mularcs) { 
      mcolor_t temp(un, arc.second, mcolor_t::Union);  
      un = temp;
    }
    return un;
  }
 
  inline bool defined(vertex_t const & v) const { 
    return (m_mularcs.count(v) != 0);
  }
 
  inline size_t size() const {
    return m_mularcs.size();  
  }

  DECLARE_ITERATOR( iter, m_mularcs, begin, begin )
  DECLARE_ITERATOR( iter, m_mularcs, end, end )
  DECLARE_CONST_ITERATOR( citer, m_mularcs, begin, cbegin )
  DECLARE_CONST_ITERATOR( citer, m_mularcs, end, cend )
  DECLARE_CONST_ITERATOR( citer, m_mularcs, cbegin, cbegin )
  DECLARE_CONST_ITERATOR( citer, m_mularcs, cend, cend )
  DECLARE_CONST_ITERATOR( criter, m_mularcs, crbegin, crbegin )
  DECLARE_CONST_ITERATOR( criter, m_mularcs, crend, crend )

private: 
  map_t m_mularcs;
}; 

}

#endif
