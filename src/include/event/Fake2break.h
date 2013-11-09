#ifndef FAKE2BREAK_H_
#define FAKE2BREAK_H_

namespace event { 

template<class mcolor_t>
struct FakeTwoBreak : public event::Event {
  typedef structure::Mularcs<mcolor_t> mularcs_t; 
  typedef std::pair<vertex_t, mcolor_t> color_arc_t;

  FakeTwoBreak(const arc_t& central, const mularcs_t& fs, const edge_t& mother)
  : fathers(fs)
  , central_edge(central)
  , mother_arc(mother) 
  {
  } 

  inline arc_t get_central_arc() const { 
    return central_edge;
  } 

  inline mcolor_t get_union_mcolor() const { 
    return mother_arc.second; 
  } 

  inline mularcs_t get_end_edges() const { 
    return fathers; 
  } 

  inline color_arc_t get_mother_edge() const { 
    return mother_arc; 
  }
private: 
  mularcs_t fathers; 
  arc_t central_edge; 
  color_arc_t mother_arc;
};

}

#endif
