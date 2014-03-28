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
  	
  TwoBreak(const arc_t& a1, const arc_t& a2, const mcolor_t& Q)
  : arcs({a1, a2}) 
  , mcolor(Q) 
  {
  }

  TwoBreak(const vertex_t& x1, const vertex_t& y1, const vertex_t& x2, const vertex_t& y2, const mcolor_t& Q)
  : arcs({arc_t(x1, y1), arc_t(x2, y2)})
  , mcolor(Q) 
  {
  }
 
  inline void change_vertex(size_t index, const vertex_t& v) {
    if (index == 0) {
      arcs[0].first = v; 
    } else if (index == 1) {
      arcs[0].second = v; 
    } else if (index == 2) {
      arcs[1].first = v; 
    } else if (index == 3) {
      arcs[1].second = v; 
    } 
  }
   
  inline TwoBreak inverse() const {
    return TwoBreak(arcs[0].first, arcs[1].first, arcs[0].second, arcs[1].second, mcolor);
  }

  inline arc_t get_arc(size_t index) const { 
    assert(index < 2); 
    return arcs[index];
  } 

  inline vertex_t get_vertex(size_t index) const { 
    if (index == 0) {
      return arcs[0].first; 
    } else if (index == 1) {
      return arcs[0].second; 
    } else if (index == 2) {
      return arcs[1].first; 
    } else if (index == 3) {
      return arcs[1].second; 
    } 
    assert(false);
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

  void apply_single(partgraph_t& local_graph) const;
private: 
  arc_t arcs[2]; // (x1,y1) x (x2,y2) = (x1,x2) + (y1,y2)
  mcolor_t mcolor; 
};

}

template<class mcolor_t>
void event::TwoBreak<mcolor_t>::apply_single(partgraph_t& local_graph) const {
  //std::cerr << arcs[0].first << " " <<  arcs[0].second << ", " << arcs[1].first << " " << arcs[1].second << std::endl; 
//" - " << genome_match::mcolor_to_name(mcolor) << std::endl;
  for(size_t i = 0; i < 2; ++i) {
    if (arcs[i].first != Infty || arcs[i].second != Infty) {
      local_graph.erase(arcs[i].first, arcs[i].second);
    }
  }
	
  if (arcs[0].first != Infty || arcs[1].first != Infty) {
    local_graph.insert(arcs[0].first, arcs[1].first);
  }

  if (arcs[0].second != Infty || arcs[1].second != Infty) {
    local_graph.insert(arcs[0].second, arcs[1].second);
  }
} 

#endif

