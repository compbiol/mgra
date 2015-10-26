//
// Created by pavel on 10/26/15.
//

#ifndef MGRA_WGD_HPP
#define MGRA_WGD_HPP

namespace event {

template<class mcolor_t>
struct WGD {

    WGD() = default;

    WGD(mcolor_t const &parent, mcolor_t const &children, size_t multiplicity)
            : m_parent(parent), m_children(children), m_multiplicity(multiplicity) {
    }

    DECLARE_GETTER(mcolor_t, m_parent, parent)

    DECLARE_GETTER(mcolor_t, m_children, children)

    DECLARE_GETTER(size_t, m_multiplicity, multiplicity)

private:
    mcolor_t m_parent;
    mcolor_t m_children;
    size_t m_multiplicity;
};

}

#endif //MGRA_WGD_HPP
