//
// Created by pavel on 10/26/15.
//

#ifndef MGRA_GENOME_GRAPH_HPP
#define MGRA_GENOME_GRAPH_HPP

#include "defined_1.hpp"

namespace graph {

template<class vertex_type>
struct GenomeGraph {
    using vertex_t = vertex_type;
    using vertex_ptr = typename vertex_t::vertex_ptr;
    using vertex_comp = typename vertex_t::vertex_comp;
    using verteces_set = std::set<vertex_ptr, vertex_comp>;

    using edge_t = std::pair<vertex_ptr, vertex_ptr>;
    using edges_set = utility::sym_multihashmap<vertex_ptr>;

    using insdel_t = event::InsDel<vertex_t, void>;
    using twobreak_t = event::TwoBreak<vertex_t, void>;

    //using clone_t =  event::Clone<structure::Mcolor>;
    //using tandem_duplication_t = event::TandemDuplication<structure::Mcolor>;

    void add_vertex(vertex_ptr v) {
        assert(vertices.count(v) == 0);
        vertices.insert(v);
    }

    void add_edge(vertex_ptr u, vertex_ptr v) {
        if (vertices.count(u) == 0) { add_vertex(u); }
        if (vertices.count(v) == 0) { add_vertex(v); }
        edges.insert(u, v);
    }

    void add_edge(edge_t const &edge) {
        add_edge(edge.first, edge.second);
    }

    void erase_vertex(vertex_ptr v) {
        assert(vertices.count(v) != 0);
        assert(!edges.defined(v));
        vertices.erase(v);
    }

    void erase_edge(vertex_ptr u, vertex_ptr v) {
        assert(vertices.count(u) != 0);
        assert(vertices.count(v) != 0);
        edges.erase(u, v);
    }

    void erase_edge(edge_t const &edge) {
        erase_edge(edge.first, edge.second);
    }

    /**
     * Apply two-break operation on genome graph. (also see event/TwoBreak.hpp)
     */
    void apply(twobreak_t const & event);

    /**
     * Apply insertion/deletion operation on genome graph. (also see event/InsDel.hpp)
     */
    void apply(insdel_t & event);

    /**
     * Apply tandem duplication operation on genome graph. (also see event/TandemDuplication.hpp)
     */
    //void apply(tandem_duplication_t & event);

    /**
     * Apply series of twobreaks on genome graph
     */
    void apply_twobreak_transformation(std::list<twobreak_t> const &transformation) {
        for (auto it = transformation.cbegin(); it != transformation.cend(); ++it) {
            apply(*it);
        }
    }

    size_t number_of_verteces() const {
        return vertices.size();
    }

    size_t number_of_edges() const {
        return edges.size();
    }

    bool is_vertex(vertex_ptr const & v) const {
        return (bool) vertices.count(v);
    }

    bool is_edge(vertex_ptr const & u, vertex_ptr const & v) const {
        assert(vertices.count(u) != 0);
        assert(vertices.count(v) != 0);
        return (edges.find(u, v) != edges.cend());
    }

    bool is_edge(edge_t const &edge) const {
        return is_edge(edge.first, edge.second);
    }

    /**
     * Count number of circular chromosomes in genome.
     * In paper about linearization algorithm this function corresponding c(.)
     * @return number of cicular chromosomes in genome.
     */
    size_t count_circular_chromosome() const;

    /**
     * Get genome from graph with <name>
     * Return genome class (see Genome.hpp)
     */
    //genome_t get_genome(std::string const &name) const;

    /**
     * Get connected components from graph.
     * Return connected components in equivalance class (see equivalence.hpp)
     */
    utility::equivalence<vertex_t> get_connected_components() const;

private:
    verteces_set vertices;
    edges_set edges;
};

template<class vertex_t>
void GenomeGraph<vertex_t>::apply(twobreak_t const &event) {
    for (size_t i = 0; i < 2; ++i) {
        if (!event.get_edge(i).first->is_infinity() || !event.get_edge(i).second->is_infinity()) {
            erase_edge(event.get_edge(i));
        }
    }

    for (size_t i = 0; i < 2; ++i) {
        if (!event.get_vertex(i)->is_infinity() || !event.get_vertex(i + 2)->is_infinity()) {
            add_edge(event.get_vertex(i), event.get_vertex(i + 2));
        }
    }
}

template <class vertex_t>
void GenomeGraph<vertex_t>::apply(insdel_t & event) {
    if (event.is_insertion()) {
        add_edge(event.get_edge());
    } else {
        erase_edge(event.get_edge());
    }
}

/*
template <class vertex_t>
void GenomeGraph<vertex_t>::apply(tandem_duplication_t & event) {
}*/

template<class vertex_t>
size_t GenomeGraph<vertex_t>::count_circular_chromosome() const {
    return 0;
}

/*template <class vertex_t>
GenomeGraph<vertex_t>::genome_t GenomeGraph<vertex_t>::get_genome(std::string const & name) const {

}

template <class vertex_t>
utility::equivalence<vertex_t> GenomeGraph<vertex_t>::get_connected_components() const {
}*/

}

#endif //MGRA_GENOME_GRAPH_HPP
