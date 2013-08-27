#ifndef MULARCS_H_
#define MULARCS_H_

namespace structure { 

template<class mcolor_t>
struct Mularcs { 
  typedef std::multimap<vertex_t, mcolor_t> map_t; 
  typedef typename map_t::const_iterator citer; 
  typedef typename map_t::const_reverse_iterator criter;
  typedef typename map_t::iterator iter;    

  inline void insert(const vertex_t& v, const mcolor_t& mc) { 
    mularcs.insert(std::make_pair(v, mc));	
  } 

  inline void insert(const vertex_t& v, size_t i) {
    if (mularcs.find(v) != mularcs.end()) { 
      mularcs.find(v)->second.insert(i);
    } else { 
      mularcs.insert(std::make_pair(v, mcolor_t(i)));
    }
  }
 
  inline void erase(const vertex_t& v) { 
    mularcs.erase(v);
  } 
 
  inline bool defined(const vertex_t& v) const { 
    return (mularcs.count(v) != 0);
  }
 
  inline vertex_t get_vertex(const mcolor_t& color) const { 
    vertex_t v = "";
    for (const auto & arc : mularcs) {
      if (arc.second == color) { 
        v = arc.first;
      }
    }
    return v;   
  } 
  
  inline mcolor_t get_multicolor(const vertex_t& v) const {
    if (mularcs.count(v) != 0) { 
      return mularcs.find(v)->second; 
    } else { 
      return mcolor_t();
    } 
  } 

  inline std::pair<citer, citer> equal_range(const vertex_t& v) const {
    return mularcs.equal_range(v);
  } 
 
  inline mcolor_t union_multicolors() const {
    mcolor_t un; 
    for (const auto &arc : mularcs) { 
      un = mcolor_t(un, arc.second, mcolor_t::Union);  
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
