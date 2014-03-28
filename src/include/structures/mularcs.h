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
    mularcs.insert(std::make_pair(v, mc));	
  } 

  inline void insert(vertex_t const & v, size_t i) {
    if (mularcs.find(v) != mularcs.end()) { 
      mularcs.find(v)->second.insert(i);
    } else { 
      mularcs.insert(std::make_pair(v, mcolor_t(i)));
    }
  }
 
  inline void erase(vertex_t const & v) { 
    mularcs.erase(v);
  } 
 
  inline bool defined(vertex_t const & v) const { 
    return (mularcs.count(v) != 0);
  }
 
  inline vertex_t get_vertex(mcolor_t const & color) const { 
    vertex_t v = "";
    for (auto const & arc : mularcs) {
      if (arc.second == color) { 
        v = arc.first;
      }
    }
    return v;   
  } 
  
  inline mcolor_t get_multicolor(vertex_t const & v) const {
    if (mularcs.count(v) != 0) { 
      return mularcs.find(v)->second; 
    } else { 
      return mcolor_t();
    } 
  } 

  inline size_t number_unique_edge() const {
    std::unordered_set<vertex_t> processed; 
    for (auto const & arc : mularcs) { 
      processed.insert(arc.first);
    }
    return processed.size();

  } 

  inline std::pair<citer, citer> equal_range(vertex_t const & v) const {
    return mularcs.equal_range(v);
  } 
 
  inline mcolor_t union_multicolors() const {
    mcolor_t un; 
    for (auto const &arc : mularcs) { 
      mcolor_t temp(un, arc.second, mcolor_t::Union);  
      un = temp;
    }
    return un;
  }

  inline std::set<mcolor_t> get_multicolors() const { 
    std::set<mcolor_t> answer; 
    for (const auto &arc : mularcs) {
      answer.insert(arc.second);
    }		
    return answer;	 
  }

  inline size_t size() const {
    return mularcs.size();  
  }

  inline iter begin() { 
    return mularcs.begin(); 
  } 

  inline iter end() { 
    return mularcs.end();
  } 

  inline citer cbegin() const { 
    return mularcs.cbegin(); 
  } 

  inline citer cend() const { 
    return mularcs.cend();
  } 

  inline criter crbegin() const { 
    return mularcs.crbegin(); 
  } 

  inline criter crend() const { 
    return mularcs.crend();
  } 
private: 
  map_t mularcs;
}; 

}
#endif
