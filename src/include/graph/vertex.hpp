//
// Created by pavel on 06.11.15.
//

#ifndef MGRA_VERTEX_HPP
#define MGRA_VERTEX_HPP

template<class info_t>
struct Vertex {
    enum tag_t {block, infinity, repeat};

    Vertex() = default;

    Vertex(size_t ind, info_t info, tag_t t)
            : index(ind), data(info), tag(t)
    {}

    Vertex(size_t ind)
            : Vertex(ind, info_t(), block)
    {}

    Vertex(size_t ind, info_t info)
            : Vertex(ind, info, block)
    {}

    bool is_infinity() const {
        return  (tag == infinity);
    }

public:
    size_t index;
    info_t data;
    tag_t tag;
};

#endif //MGRA_VERTEX_HPP
