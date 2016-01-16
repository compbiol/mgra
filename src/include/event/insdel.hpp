//
// Created by pavel on 1/15/16.
//

#ifndef MGRA_INSDEL_HPP
#define MGRA_INSDEL_HPP

namespace event {

template<class vertex_type, class mcolor_t>
struct InsDel : public InsDel<vertex_type, void> {
    using citer = typename mcolor_t::citer;

    using basic_insdel_t = InsDel<vertex_type, void>;
    using edge_t = typename basic_insdel_t::edge_t;
    using vertex_t = typename basic_insdel_t::vertex_t;
    using vertex_ptr = typename basic_insdel_t::vertex_ptr;
    using dependence_type = typename basic_insdel_t::dependence_type;

    using insdel_t = InsDel<vertex_t, mcolor_t>;
    using insdel_ptr = std::shared_ptr<InsDel<vertex_t, mcolor_t>>;

    InsDel(vertex_t const &x, vertex_t const &y, mcolor_t const &multicolor, bool is_insertion = true)
    : basic_insdel_t(x, y, is_insertion)
    , m_multicolor(multicolor)
    {
        assert(!m_multicolor.empty());
    }

    InsDel(edge_t const &edge, mcolor_t const &multicolor, bool is_insertion = true)
    : InsDel(edge.first, edge.second, multicolor, is_insertion)
    { }

    InsDel(basic_insdel_t const & insdel, mcolor_t const & multicolor)
    : InsDel(insdel.get_edge().first, insdel.get_edge().second, multicolor)
    { }

    DECLARE_GETTER(mcolor_t, m_multicolor, multicolor)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, begin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, end, cend)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, cbegin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, cend, cend)

protected:
    mcolor_t m_multicolor;
};

template<class vertex_type>
struct InsDel<vertex_type, void> {
    using insdel_t = InsDel<vertex_type, void>;
    using insdel_ptr = std::shared_ptr<InsDel>;
    using vertex_t = vertex_type;
    using vertex_ptr = typename vertex_t::vertex_ptr;
    using edge_t = std::pair<vertex_ptr, vertex_ptr>;

    InsDel(vertex_ptr const &x, vertex_ptr const &y, bool is_insertion = true)
    : m_edge(x, y)
    , m_is_insertion(is_insertion)
    {
        assert(!x->is_infinity());
        assert(!y->is_infinity());
    }

    InsDel(edge_t const &edge, bool is_insertion = true)
    : InsDel(edge.first, edge.second, is_insertion)
    { }

    template<class mcolor_t>
    explicit InsDel(InsDel<vertex_t, mcolor_t> const & insdel)
    : InsDel(insdel.m_edge.first, insdel.m_edge.second, insdel.m_is_insertion)
    { }

    bool is_insertion() const {
        return m_is_insertion;
    }

    DECLARE_GETTER(edge_t, m_edge, edge)

protected:
    edge_t m_edge; // (x1, y1)
    bool m_is_insertion;
};

}

#endif //MGRA_INSDEL_HPP
