#ifndef MULARCS_H_
#define MULARCS_H_

template<class mcolor_t>
struct Mularcs { 
  typedef std::multimap<vertex_t, mcolor_t> mymap; 
  typedef typename mymap::const_iterator citer; 
  typedef typename mymap::const_reverse_iterator criter;
  typedef typename mymap::iterator iter;    

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
 
  inline mcolor_t union_multicolors() const {
    mcolor_t un; 
    for (auto it = mularcs.cbegin(); it != mularcs.cend(); ++it) {
      un = mcolor_t(un, it->second, mcolor_t::Union);  
    }
    return un;
  }
 
  inline mcolor_t get_multicolor(const vertex_t& v) const {
    if (mularcs.find(v) != mularcs.cend()) { 
      return mularcs.find(v)->second; 
    } else { 
      return mcolor_t();
    } 
  } 

  inline std::set<mcolor_t> get_multicolors() const { 
    std::set<mcolor_t> answer; 
    for (auto it = mularcs.cbegin(); it != mularcs.cend(); ++it) {
      answer.insert(it->second);
    }		
    return answer;	 
  }

  inline iter find(const vertex_t& v) { 
    return mularcs.find(v);
  }

  inline citer find(const vertex_t& v) const { 
    return mularcs.find(v);
  }
 
  inline size_t size() const {
    return mularcs.size();  
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
  mymap mularcs;
}; 

#endif
