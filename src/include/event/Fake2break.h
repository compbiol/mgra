#ifndef FAKE2BREAK_H_
#define FAKE2BREAK_H_

namespace event { 

template<class mcolor_t>
struct FakeTwoBreak : public event::Event {
  typedef structure::Mularcs<mcolor_t> mularcs_t; 
  typedef std::pair<vertex_t, mcolor_t> color_arc_t;

  FakeTwoBreak() 
  : m_with_pseudo_vertex(false)
  {
  }

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

  inline void change_central_edge_first(vertex_t const & v) {
    central_edge.first = v; 
  }

  inline void change_central_edge_second(vertex_t const & v) {
    central_edge.second = v; 
  }

  inline void change_mother_vertex(vertex_t const & v) {
    mother_arc.first = v;
  }

  DECLARE_GETTER( vertex_t, central_arc.first, father_vertex )
  DECLARE_GETTER( arc_t, central_edge, central_arc )
  DECLARE_GETTER( mularcs_t, fathers, end_edges )
  DECLARE_GETTER( color_arc_t, mother_arc, mother_edge )
  DECLARE_GETTER( mcolor_t, mother_arc.second, mcolor )

private: 
  mularcs_t fathers; 
  arc_t central_edge; 
  color_arc_t mother_arc;
  bool m_with_pseudo_vertex;
};

}

#endif
