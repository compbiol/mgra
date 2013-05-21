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
 
  inline bool is_simple_vertice() const {
    if (mularcs.size() == 2 && mularcs.cbegin()->second.is_good_multiedge() && mularcs.crbegin()->second.is_good_multiedge()) { 
      return true; 
    } 
    return false; 
  }  
 
  inline bool is_duplication_vertice() const {	
    for(auto im = mularcs.cbegin(); im != mularcs.cend(); ++im) { 
      for(auto it = mularcs.cbegin(); it != mularcs.cend(); ++it) {
	if (*im == *it) { 
	  continue;
	} 
	mcolor_t color(im->second, it->second, mcolor_t::Intersection);
	if (!color.empty()) { 
	  return true; 
	} 
      } 
    }  
    return false; 
  } 

  inline bool is_fair_vertice() const {
    if (mularcs.size() == 3) {
      for(auto it = mularcs.cbegin(); it != mularcs.cend(); ++it) {
	if (!it->second.is_good_multiedge()) { 
	  return false;
	}  
      } 
      return true; 
    }  
    return false; 
  }

  inline mcolor_t get_multicolor(const vertex_t& v) const {
    if (mularcs.find(v) != mularcs.cend()) { 
      return mularcs.find(v)->second; 
    } else { 
      return mcolor_t();
    } 
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
