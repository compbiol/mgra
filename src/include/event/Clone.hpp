#ifndef MGRA_CLONE_HPP
#define MGRA_CLONE_HPP

namespace event {

template<class mcolor_t>
struct Clone {
    using mularcs_t = structure::Mularcs<mcolor_t>;
    using arc_t = std::pair<vertex_t, mcolor_t>;
    using edge_t = std::pair<vertex_t, vertex_t>;

    Clone()
            : m_with_pseudo_vertex(false) {
    }

    Clone(edge_t const &central_edge, mularcs_t const &fathers, arc_t const &mother_arc, bool with_pseudo_vertex)
            : m_fathers(fathers), m_central_edge(central_edge), m_mother_arc(mother_arc),
              m_with_pseudo_vertex(with_pseudo_vertex) {
    }

    inline bool is_have_pseudo_vertex() const {
        return m_with_pseudo_vertex;
    }

    inline arc_t get_mother_arc() const {
        return m_mother_arc;
    }

    inline edge_t get_central_edge() const {
        return m_central_edge;
    }

    inline mularcs_t get_fathers() const {
        return m_fathers;
    }

    DECLARE_GETTER(vertex_t, central_arc.first, father_vertex)

    DECLARE_GETTER(edge_t, m_central_edge, central_arc)

    DECLARE_GETTER(mularcs_t, m_fathers, end_edges)

    DECLARE_GETTER(arc_t, m_mother_arc, mother_edge)

    DECLARE_GETTER(mcolor_t, m_mother_arc.second, mcolor)

private:
    mularcs_t m_fathers;
    edge_t m_central_edge;
    arc_t m_mother_arc;
    bool m_with_pseudo_vertex;
};

}

#endif // MGRA_CLONE_HPP
