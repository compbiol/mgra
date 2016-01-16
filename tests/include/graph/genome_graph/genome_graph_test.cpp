//
// Created by pavel on 12/29/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"
#include "graph/vertex.hpp"
#include "graph/genome_graph/genome_graph.hpp"

using namespace std;
using graph_t = graph::GenomeGraph<Vertex<std::string>>;

using edge_t = graph_t::edge_t;
using vertex_t = graph_t::vertex_t;
using vertex_ptr = graph_t::vertex_ptr;
using twobreak_t = graph_t::twobreak_t;

class GenomeGraphTest : public testing::Test {

    virtual void SetUp() override {
        verteces.push_back(vertex_ptr(new vertex_t(1, "1h"))); //0
        verteces.push_back(vertex_ptr(new vertex_t(2, "1t"))); //1
        verteces.push_back(vertex_ptr(new vertex_t(3, "2h"))); //2
        verteces.push_back(vertex_ptr(new vertex_t(4, "2t"))); //3
        verteces.push_back(vertex_ptr(new vertex_t(5, "3t"))); //4
        verteces.push_back(vertex_ptr(new vertex_t(6, "3h"))); //5
        verteces.push_back(vertex_ptr(new vertex_t(7, "3t"))); //6
        verteces.push_back(vertex_ptr(new vertex_t(8, "3h"))); //7

        //add zero vertex from add_edge function
        graph.add_vertex(verteces[6]);
        graph.add_vertex(verteces[7]);
        graph.add_edge(edge_t(verteces[6], verteces[7]));
        graph.add_edge(verteces[6], verteces[7]);

        //add one vertex from add_edge function
        graph.add_vertex(verteces[0]);
        graph.add_vertex(verteces[1]);
        graph.add_edge(edge_t(verteces[1], verteces[3]));
        graph.add_edge(verteces[0], verteces[2]);

        //add two vertex from add_edge function
        graph.add_edge(edge_t(verteces[4], verteces[5]));
        graph.add_edge(verteces[4], verteces[5]);
    }

protected:
    std::vector<vertex_ptr> verteces;
    graph_t graph;
};

TEST_F(GenomeGraphTest, AddFunction) {
    EXPECT_EQ(graph.number_of_verteces(), verteces.size());
    for (size_t i = 0; i < verteces.size(); ++i) {
        EXPECT_TRUE(graph.is_vertex(verteces[i]));
    }

    EXPECT_EQ(graph.number_of_edges(), 6);
    EXPECT_TRUE(graph.is_edge(verteces[6], verteces[7]));
    EXPECT_TRUE(graph.is_edge(verteces[1], verteces[3]));
    EXPECT_TRUE(graph.is_edge(edge_t(verteces[0], verteces[2])));
    EXPECT_TRUE(graph.is_edge(edge_t(verteces[4], verteces[5])));
}

TEST_F(GenomeGraphTest, EraseFunction) {
    graph.erase_edge(verteces[6], verteces[7]);
    EXPECT_EQ(graph.number_of_edges(), 5);
    EXPECT_TRUE(graph.is_edge(verteces[6], verteces[7]));

    graph.erase_edge(edge_t(verteces[6], verteces[7]));
    EXPECT_EQ(graph.number_of_edges(), 4);
    EXPECT_FALSE(graph.is_edge(edge_t(verteces[6], verteces[7])));

    graph.erase_edge(edge_t(verteces[1], verteces[3]));
    EXPECT_EQ(graph.number_of_edges(), 3);
    EXPECT_FALSE(graph.is_edge(verteces[1], verteces[3]));

    graph.erase_edge(verteces[0], verteces[2]);
    EXPECT_EQ(graph.number_of_edges(), 2);
    EXPECT_FALSE(graph.is_edge(edge_t(verteces[0], verteces[2])));

    graph.erase_edge(verteces[4], verteces[5]);
    EXPECT_EQ(graph.number_of_edges(), 1);
    EXPECT_TRUE(graph.is_edge(edge_t(verteces[4], verteces[5])));

    graph.erase_edge(edge_t(verteces[4], verteces[5]));
    EXPECT_EQ(graph.number_of_edges(), 0);
    EXPECT_FALSE(graph.is_edge(verteces[4], verteces[5]));

    for (size_t i = 0; i < verteces.size(); ++i) {
        graph.erase_vertex(verteces[i]);
        EXPECT_FALSE(graph.is_vertex(verteces[i]));
        EXPECT_EQ(graph.number_of_verteces(), verteces.size() - i - 1);
    }
}

TEST_F(GenomeGraphTest, ApplyTwoBreak) {
    graph.apply(twobreak_t(verteces[6], verteces[7], verteces[4], verteces[5]));
    EXPECT_EQ(graph.number_of_edges(), 6);
    EXPECT_EQ(graph.number_of_verteces(), verteces.size());
    EXPECT_TRUE(graph.is_edge(verteces[6], verteces[4]));
    EXPECT_TRUE(graph.is_edge(verteces[5], verteces[7]));

    //vertex_ptr v(new vertex_t(9, "oo", vertex_t::infinity));
    //TODO with infinity start to work
}