//
// Created by pavel on 11/6/15.
//

#ifndef MGRA_TWO_BREAK_HPP
#define MGRA_TWO_BREAK_HPP

namespace event {

template<class vertex_t, class mcolor_t>
struct TwoBreak {
    enum dependence_type {
        independent, weakly_dependent, strong_dependent
    };

    using citer = typename mcolor_t::citer;
    using vertex_ptr = typename vertex_t::vertex_ptr;
    using edge_t = std::pair<vertex_ptr, vertex_ptr>;
    using twobreak_ptr = std::shared_ptr <TwoBreak<vertex_t, mcolor_t>>;

    TwoBreak() = default;

    TwoBreak(vertex_ptr const &x1, vertex_ptr const &y1, vertex_ptr const &x2, vertex_ptr const &y2,
             mcolor_t const &multicolor)
            : m_multicolor(multicolor) {
        m_verteces[0] = x1;
        m_verteces[1] = y1;
        m_verteces[2] = x2;
        m_verteces[3] = y2;
    }

    TwoBreak(edge_t const &a1, edge_t const &a2, mcolor_t const &multicolor)
            : TwoBreak(a1.first, a1.second, a2.first, a2.second, multicolor) { }

    inline void change_vertex(size_t index, vertex_ptr const &v) {
        assert(index < 4);
        m_verteces[index] = v;
    }

    inline TwoBreak inverse() const {
        return TwoBreak(m_verteces[0], m_verteces[2], m_verteces[1], m_verteces[3], m_multicolor);
    }

    inline edge_t get_edge(size_t index) const {
        assert(index < 2);
        if (index == 0) {
            return edge_t(m_verteces[0], m_verteces[1]);
        }
        return edge_t(m_verteces[2], m_verteces[3]);
    }

    inline vertex_ptr get_vertex(size_t index) const {
        assert(index < 4);
        return m_verteces[index];
    }


    inline mcolor_t get_multicolor() const {
        return m_multicolor;
    }

    /*
     * This function tested two two-break on depeneds.
     * Return: indepented, weakly dependent, strong dependent
     */
    dependence_type is_dependent(TwoBreak const &tested) const;

    DECLARE_GETTER(mcolor_t, m_multicolor, mcolor)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, begin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, end, cend)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, cbegin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_multicolor, cend, cend)

private:
    std::array<vertex_ptr, 4> m_verteces; // (x1,y1) x (x2,y2) = (x1,x2) x (y1,y2)
    mcolor_t m_multicolor;
};

template<class vertex_t>
struct TwoBreak<vertex_t, void> {
    enum dependence_type {
        independent, weakly_dependent, strong_dependent
    };

    using twobreak_ptr = std::shared_ptr<TwoBreak>;
    using vertex_ptr = typename vertex_t::vertex_ptr;
    using edge_t = std::pair<vertex_ptr, vertex_ptr>;

    TwoBreak() = default;

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

    inline void change_vertex(size_t index, vertex_ptr const &v) {
        assert(index < 4);
        m_verteces[index] = v;
    }

    inline TwoBreak inverse() const {
        return TwoBreak(m_verteces[0], m_verteces[2], m_verteces[1], m_verteces[3]);
    }

    bool operator==(TwoBreak const & b) const {
        return (m_verteces == b.m_verteces);
    }

    inline edge_t get_edge(size_t index) const {
        assert(index < 2);
        if (index == 0) {
            return edge_t(m_verteces[0], m_verteces[1]);
        }
        return edge_t(m_verteces[2], m_verteces[3]);
    }

    inline vertex_ptr get_vertex(size_t index) const {
        assert(index < 4);
        return m_verteces[index];
    }

    /*
     * This function tested two two-break on depeneds.
     * Return: indepented, weakly dependent, strong dependent
     */
    dependence_type is_dependent(TwoBreak const &tested) const;

private:
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

template<class vertex_t, class mcolor_t>
typename TwoBreak<vertex_t, mcolor_t>::dependence_type TwoBreak<vertex_t, mcolor_t>::is_dependent(
        TwoBreak<vertex_t, mcolor_t> const &tested) const {
    mcolor_t inter_color(m_multicolor, tested.m_multicolor, mcolor_t::Intersection);

    if (inter_color.empty()) {
        return independent;
    } else {
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

}
#endif //MGRA_TWO_BREAK_HPP
