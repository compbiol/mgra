//
// Created by pavel on 10/26/15.
//

#ifndef MGRA_GENOME_GRAPH_HPP
#define MGRA_GENOME_GRAPH_HPP

#define "defined.hpp"

template<class vertex_t>
struct GenomeGraph {
    using genome_t = structure::Genome;
    using twobreak_t = event::TwoBreak<structure::Mcolor>;
    using clone_t =  event::Clone<structure::Mcolor>;
    using insdel_t = event::InsDel<structure::Mcolor>;
    using tandem_duplication_t = event::TandemDuplication<structure::Mcolor>;


    inline void add_vertex(vertex_t const & v) {
        assert(vertices.count(v) != 0);
        vertices.insert(v);
    }

    inline void add_edge(vertex_t const & u, vertex_t const & v) {
        if (vertices.count(u) != 0) {add_vertex(u);}
        if (vertices.count(v) != 0) {add_vertex(v);}
        edges.insert(u, v);
    }

    inline void erase_vertex(vertex_t const & v) {
        assert(vertices.count(v) == 0);
        vertices.erase(v);
    }

    inline void erase_edge(vertex_t const & u, vertex_t const & v) {
        assert(vertices.count(u) == 0);
        assert(vertices.count(v) == 0);
        edges.erase(u, v);
    }

    /**
     * Apply two-break operation on genome graph. (also see event/TwoBreak.hpp)
     */
    void apply(twobreak_t & event);

    /**
     * Apply clone operation on genome graph. (also see event/Clone.hpp. About clone operations see in mgra2 paper)
     */
    void apply(clone_t & event);

    /**
     * Apply insertion/deletion operation on genome graph. (also see event/InsDel.hpp)
     */
    void apply(insdel_t & event);

    /**
     * Apply tandem duplication operation on genome graph. (also see event/TandemDuplication.hpp)
     */
    void apply(tandem_duplication_t & event);

    /**
     * Apply series of twobreaks on genome graph
     */
    void apply_twobreaks_transformation(std::list<twobreak_t> const & transformation) {
        for (auto it = transformation.cbegin(); it != transformation.cend(); ++it) {
            apply(*it);
        }
    }

    /**
     * Count number of circular chromosomes in genome.
     * In paper about linearization algorithm this function corresponding c(.)
     * Return number of cicular chromosomes in genome.
     */
    size_t count_circular_chromosome() const;

    /**
     * Get genome from graph with <name>
     * Return genome class (see Genome.hpp)
     */
    genome_t get_genome(std::string const & name) const;

    /**
     * Get connected components from graph.
     * Return connected components in equivalance class (see equivalence.hpp)
     */
    utility::equivalence<vertex_t> get_connected_components() const;

private:
    std::unordered_set<vertex_t> vertices;
    utility::sym_multihashmap<vertex_t> edges;
};

template <class vertex_t>
void GenomeGraph<vertex_t>::apply(twobreak_t & event) {
}

template <class vertex_t>
void GenomeGraph<vertex_t>::apply(clone_t & event) {
}

template <class vertex_t>
void GenomeGraph<vertex_t>::apply(insdel_t & event) {
}

template <class vertex_t>
void GenomeGraph<vertex_t>::apply(tandem_duplication_t & event) {
}

template <class vertex_t>
size_t GenomeGraph<vertex_t>::count_circular_chromosome() const {
}

template <class vertex_t>
GenomeGraph<vertex_t>::genome_t GenomeGraph<vertex_t>::get_genome(std::string const & name) const {
}

template <class vertex_t>
utility::equivalence<vertex_t> GenomeGraph<vertex_t>::get_connected_components() const {
}

#endif //MGRA_GENOME_GRAPH_HPP
