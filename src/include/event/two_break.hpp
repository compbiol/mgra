//
// Created by pavel on 11/6/15.
//

#ifndef MGRA_TWO_BREAK_HPP
#define MGRA_TWO_BREAK_HPP

namespace event {

template<class vertex_type, class mcolor_t>
struct TwoBreak : public TwoBreak<vertex_type, void> {
    using basic_twobreak_t = TwoBreak<vertex_type, void>;
    using edge_t = typename basic_twobreak_t::edge_t;
    using vertex_t = typename basic_twobreak_t::vertex_t;
    using vertex_ptr = typename basic_twobreak_t::vertex_ptr;
    using dependence_type = typename basic_twobreak_t::dependence_type;
    using citer = typename mcolor_t::citer;
    using twobreak_t = TwoBreak<vertex_t, mcolor_t>;
    using twobreak_ptr = std::shared_ptr<TwoBreak<vertex_t, mcolor_t>>;

    TwoBreak(vertex_ptr x1, vertex_ptr y1, vertex_ptr x2, vertex_ptr y2, mcolor_t const &multicolor)
    : basic_twobreak_t(x1, y1, x2, y2)
    , m_multicolor(multicolor)
    {
        assert(!m_multicolor.empty());
    }

    TwoBreak(edge_t const &a1, edge_t const &a2, mcolor_t const &multicolor)
    : basic_twobreak_t(a1.first, a1.second, a2.first, a2.second)
    , m_multicolor(multicolor)
    {
        assert(!m_multicolor.empty());
    }

    TwoBreak(basic_twobreak_t const & twobreak, mcolor_t const &multicolor)
    : basic_twobreak_t(twobreak)
    , m_multicolor(multicolor)
    {
        assert(!m_multicolor.empty());
    }

    TwoBreak inverse() const {
        return TwoBreak(basic_twobreak_t::inverse(), m_multicolor);
    }

    mcolor_t get_multicolor() const {
        return m_multicolor;
    }

    /*
     * This function tested two two-break on depeneds.
     * Return: indepented, weakly dependent, strong dependent
     */
    dependence_type is_dependent(TwoBreak const &tested) const {
        mcolor_t inter_color(m_multicolor, tested.m_multicolor, mcolor_t::Intersection);
        if (inter_color.empty()) {
            return dependence_type::independent;
        } else {
            return basic_twobreak_t::is_dependent(tested);
        }
    }

    DECLARE_GETTER(mcolor_t, m_multicolor, mcolor)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, begin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, end, cend)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, cbegin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, cend, cend)

protected:
    mcolor_t m_multicolor;
};

template<class vertex_type>
struct TwoBreak<vertex_type, void> {
    enum dependence_type {
        independent, weakly_dependent, strong_dependent
    };

    using twobreak_t = TwoBreak<vertex_type, void>;
    using twobreak_ptr = std::shared_ptr<TwoBreak>;
    using vertex_t = vertex_type;
    using vertex_ptr = typename vertex_t::vertex_ptr;
    using edge_t = std::pair<vertex_ptr, vertex_ptr>;

    TwoBreak(vertex_ptr const &x1, vertex_ptr const &y1, vertex_ptr const &x2, vertex_ptr const &y2) {
        m_verteces[0] = x1;
        m_verteces[1] = y1;
        m_verteces[2] = x2;
        m_verteces[3] = y2;
    }

    TwoBreak(edge_t const &a1, edge_t const &a2)
    : TwoBreak(a1.first, a1.second, a2.first, a2.second)
    { }

    template<class mcolor_t>
    explicit TwoBreak(TwoBreak<vertex_t, mcolor_t> const &two_break)
    : TwoBreak(two_break.get_vertex(0), two_break.get_vertex(1), two_break.get_vertex(2), two_break.get_vertex(3))
    { }

    void change_vertex(size_t index, vertex_ptr const &v) {
        assert(index < 4);
        m_verteces[index] = v;
    }

    twobreak_t inverse() const {
        return twobreak_t(m_verteces[0], m_verteces[2], m_verteces[1], m_verteces[3]);
    }

    bool operator==(twobreak_t const & b) const {
        return (m_verteces == b.m_verteces);
    }

    edge_t get_edge(size_t index) const {
        assert(index < 2);
        if (index == 0) {
            return edge_t(m_verteces[0], m_verteces[1]);
        }
        return edge_t(m_verteces[2], m_verteces[3]);
    }

    vertex_ptr get_vertex(size_t index) const {
        assert(index < 4);
        return m_verteces[index];
    }

    /*
     * This function tested two two-break on depeneds.
     * Return: indepented, weakly dependent, strong dependent
     */
    virtual dependence_type is_dependent(TwoBreak const &tested) const;

protected:
    std::array<vertex_ptr, 4> m_verteces; // (x1,y1) x (x2,y2) = (x1,x2) x (y1,y2)
};

template<class vertex_t>
typename TwoBreak<vertex_t, void>::dependence_type TwoBreak<vertex_t, void>::is_dependent(
        TwoBreak<vertex_t, void> const &tested) const {
    auto check_lambda = [&](size_t ind1, size_t ind2, size_t ind3) -> size_t {
        if (!(tested.m_verteces[ind1]->is_infinity() && tested.m_verteces[ind2]->is_infinity())) {
            if (tested.get_edge(ind1) == edge_t(m_verteces[ind2], m_verteces[ind3])) {
                return 1;
            }
            if (tested.get_edge(ind1) == edge_t(m_verteces[ind3], m_verteces[ind2])) {
                return 1;
            }
        }
        return 0;
    };

    size_t first = check_lambda(0, 0, 2) + check_lambda(1, 1, 3);
    size_t second = check_lambda(0, 1, 3) + check_lambda(1, 0, 2);

    if (std::max(first, second) == 1) {
        return weakly_dependent;
    } else if (std::max(first, second) == 2) {
        return strong_dependent;
    } else {
        return independent;
    }
}

}
#endif //MGRA_TWO_BREAK_HPP
