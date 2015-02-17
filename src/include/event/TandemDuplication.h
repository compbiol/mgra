#ifndef TANDEMDUPLICATION_H_
#define TANDEMDUPLICATION_H_

namespace event { 

template <class mcolor_t>
struct TandemDuplication {
  using citer = typename mcolor_t::citer;
  using edge_t = std::pair<vertex_t, vertex_t>;
  using citer_vect = std::vector<edge_t>::const_iterator;  
	
  TandemDuplication(const std::vector<edge_t>& es, const mcolor_t& Q, bool is_del, bool is_revs)
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
  std::vector<edge_t> edges; 
  mcolor_t mcolor; 
  bool is_reverse;
  bool is_deletion; 
}; 

} 

#endif 
