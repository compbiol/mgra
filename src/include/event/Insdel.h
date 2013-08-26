#ifndef INSDEL_H_
#define INSDEL_H_

namespace event { 

template<class mcolor_t>
struct InsDel {
  typedef std::pair<vertex_t, vertex_t> arc_t;
  typedef typename mcolor_t::citer citer; 
	
  InsDel(const arc_t& e, const mcolor_t& Q, bool is_del)
  : edge(e) 
  , mcolor(Q) 
  , is_deletion(is_del)
  {
  }

  InsDel(const vertex_t& x1, const vertex_t& y1, const mcolor_t& Q, bool is_del)
  : edge(arc_t(x1, y1)) 
  , mcolor(Q) 
  , is_deletion(is_del)
  {
  }
    
  inline InsDel inverse() const {
    return InsDel(edge, mcolor, !is_deletion);
  }

  inline bool is_deletion_oper() const { 
    return is_deletion;
  }
 
  inline arc_t get_edge() const { 
    return edge;
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

private: 
  arc_t edge; // (x1, y1)
  mcolor_t mcolor; 
  bool is_deletion; //FIXME: remove, only insertion
};

}
#endif

