//
// Created by pavel on 06.11.15.
//

#ifndef MGRA_VERTEX_HPP
#define MGRA_VERTEX_HPP

template<class vertex_t>
struct VertexPtrComporator {
    using vertex_ptr = typename vertex_t::vertex_ptr;

    bool operator() (vertex_ptr const &v1, vertex_ptr const &v2) const {
        return (v1->index < v2->index);
    }
};

template<class vertex_t>
struct VertexPtrHash {
    using vertex_ptr = typename vertex_t::vertex_ptr;

    std::size_t operator()(vertex_ptr const& v) const {
        return std::hash<size_t>()(v->index);
    }
};

template<class info_t>
struct Vertex {
    using vertex_ptr = std::shared_ptr<Vertex>;
    using vertex_comp = VertexPtrComporator<Vertex>;
    using vertex_hash = VertexPtrHash<Vertex>;

    enum tag_t {
        block, infinity, repeat
    };

    Vertex() = delete;

    explicit Vertex(size_t ind)
    : Vertex(ind, info_t(), block)
    { }

    Vertex(size_t ind, info_t info)
    : Vertex(ind, info, block)
    { }

    Vertex(size_t ind, info_t info, tag_t t)
    : index(ind), data(info), tag(t)
    {
        assert(ind != 0);
    }

    bool is_infinity() const {
        return  (tag == infinity);
    }

    bool operator==(Vertex v) const {
        return (index == v.index);
    }

public:
    size_t index;
    info_t data;
    tag_t tag;
};


#endif //MGRA_VERTEX_HPP
