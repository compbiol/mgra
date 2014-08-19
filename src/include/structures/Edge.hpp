#ifndef EDGE_TYPE_HPP
#define EDGE_TYPE_HPP

template<class vertex_t>
struct Edge { 
  Edge() 
  : m_is_first_mother(false)
  , m_is_second_mother(false)
  {
  }

  Edge(vertex_t const & v, mcolor_t const & Q, vertex_t const & u)
  : first(v)
  , color(Q)
  , second(u)
  , m_is_first_mother(false) 
  , m_is_second_mother(false)
  {
  }

  bool is_real() const { 
    return (!m_is_first_mother && !m_is_second_mother);
  }

public: 
  vertex_t first; 
  mcolor_t color; 
  vertex_t second;

private: 
  bool m_is_first_mother;
  mcolor_t clone_first_color; 

  bool m_is_second_mother;
  mcolor_t clone_second_color;   
};

#endif 