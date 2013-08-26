#ifndef TANDEMDUPLICATION_H_
#define TANDEMDUPLICATION_H_

namespace event { 

template <class mcolor_t>
struct TandemDuplication {
  typedef std::pair<vertex_t, vertex_t> arc_t;
  typedef typename mcolor_t::citer citer;
  typedef std::vector<arc_t>::const_iterator citer_vect;  
	
  TandemDuplication(const std::vector<arc_t>& es, const mcolor_t& Q, bool is_del, bool is_revs)
  : edges(es) 
  , mcolor(Q) 
  , is_reverse(is_revs)
  , is_deletion(is_del)
  {
  }
    
  inline TandemDuplication inverse() const {
    return TandemDuplication(edges, mcolor, !is_deletion);
  }

  inline bool is_reverse_tandem_duplication() const {
    return is_reverse;
  }

  inline bool is_deletion_oper() const { 
    return is_deletion;
  }
 
  inline mcolor_t get_mcolor() const { 
    return mcolor; 
  } 

  inline citer begin() const { 
    return mcolor.cbegin();
  } 

  inline citer end() const { 
    return mcolor.cend();
  } 

  inline citer_vect cbegin_edges() const {
    return edges.cbegin();
  }

  inline citer_vect cend_edges() const {
    return edges.cend();
  }

private: 
  std::vector<arc_t> edges; 
  mcolor_t mcolor; 
  bool is_reverse;
  bool is_deletion; 
}; 

} 

#endif 
