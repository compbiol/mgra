//
// Created by pavel on 10/26/15.
//

#ifndef MGRA_WGD_HPP
#define MGRA_WGD_HPP

namespace event {

template<class info_t>
struct wgd {
    using citer = typename std::vector<info_t>::const_iterator;
    using criter = typename std::vector<info_t>::const_reverse_iterator;


    wgd() = delete;

    wgd(info_t const &parent, info_t const &children, std::vector<info_t> const & names)
    : m_parent(parent), m_children(children), m_names(names) {
        assert(!parent.empty());
        assert(!children.empty());
        assert(names.size() >= 1);
    }

    DECLARE_GETTER(info_t, m_parent, parent)

    DECLARE_GETTER(info_t, m_children, children)

    DECLARE_CONST_ITERATOR(citer, m_names, begin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_names, end, cend)

    DECLARE_CONST_ITERATOR(citer, m_names, cbegin, cbegin)

    DECLARE_CONST_ITERATOR(citer, m_names, cend, cend)

    DECLARE_CONST_ITERATOR(criter, m_names, crbegin, crbegin)

    DECLARE_CONST_ITERATOR(criter, m_names, crend, crend)
private:
    info_t m_parent;
    info_t m_children;
    std::vector<info_t> m_names;
};

}

#endif //MGRA_WGD_HPP
