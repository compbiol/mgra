#ifndef MGRA_INSDEL_HPP
#define MGRA_INSDEL_HPP

namespace event {

template<class mcolor_t>
struct InsDel {
    using citer = typename mcolor_t::citer;
    using edge_t = std::pair<vertex_t, vertex_t>;

    InsDel(edge_t const &edge, mcolor_t const &multicolor, bool is_insertion = true)
            : m_edge(edge), m_multicolor(multicolor), m_is_insertion(is_insertion) {
    }

    InsDel(vertex_t const &x, vertex_t const &y, mcolor_t const &multicolor, bool is_insertion = true)
            : m_edge(edge_t(x, y)), m_multicolor(multicolor), m_is_insertion(is_insertion) {
    }

    bool is_insertion() const {
        return m_is_insertion;
    }

    DECLARE_GETTER(edge_t, m_edge, edge)

    DECLARE_GETTER(mcolor_t, m_multicolor, mcolor)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, begin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, end, cend)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, cbegin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, cend, cend)

private:
    edge_t m_edge; // (x1, y1)
    mcolor_t m_multicolor;
    bool m_is_insertion;
};

}

#endif //MGRA_INSDEL_HPP

