//
// Created by pavel on 10/26/15.
//

#ifndef MGRA_GENOME_GRAPH_HPP
#define MGRA_GENOME_GRAPH_HPP

//TODO
struct GenomeGraph {

    void add_vertex(vertex_t const & v);

    void add_edge(vertex_t const & u, vertex_t const & v);

    void erase_vertex(vertex_t const & v);

    void erase_edge(vertex_t const & u, vertex_t const & v);

    //apply(twobreak_t & event)
    //apply(clone_t & event)
    //apply()
    //apply()

    //get_components

    size_t count_circular_chromosome() const;
};

#endif //MGRA_GENOME_GRAPH_HPP
