#ifndef FAKE2BREAK_H_
#define FAKE2BREAK_H_

namespace event { 

template<class mcolor_t>
struct FakeTwoBreak : public event::Event {
  typedef structure::Mularcs<mcolor_t> mularcs_t; 
  typedef std::pair<vertex_t, mcolor_t> color_arc_t;

  FakeTwoBreak(arc_t const & central, mularcs_t const & fs, edge_t const & mother, bool with_pseudo_vertex)
  : fathers(fs)
  , central_edge(central)
  , mother_arc(mother) 
  , m_with_pseudo_vertex(with_pseudo_vertex)
  {
  } 

  inline bool is_have_pseudo_vertex() const {
    return m_with_pseudo_vertex;
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
  bool m_with_pseudo_vertex;
};

}

#endif
