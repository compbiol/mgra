//
// Created by pavel on 11/6/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"
#include "graph/vertex.hpp"

using namespace std;

using vertex_t = Vertex<std::string>;
using vertex_ptr = vertex_t::vertex_ptr;

TEST(VertexTest, ConstructorWithOneParameter) {
    vertex_t v(1);
    EXPECT_EQ(1, v.index);
    EXPECT_EQ(std::string(), v.data);
    EXPECT_EQ(vertex_t::tag_t::block, v.tag);
    EXPECT_TRUE(!v.is_infinity());
    EXPECT_EQ(v, v);
}

TEST(VertexTest, ConstructorWithTwoParameter) {
    vertex_t v(1, "1h");
    EXPECT_EQ(1, v.index);
    EXPECT_EQ(std::string("1h"), v.data);
    EXPECT_EQ(vertex_t::tag_t::block, v.tag);
    EXPECT_TRUE(!v.is_infinity());
}

TEST(VertexTest, ConstructorWithThreeParameter) {
    vertex_t v(1, "1h", vertex_t::tag_t::block);
    EXPECT_EQ(1, v.index);
    EXPECT_EQ(std::string("1h"), v.data);
    EXPECT_EQ(vertex_t::tag_t::block, v.tag);
    EXPECT_TRUE(!v.is_infinity());

    vertex_t v1(2, "1h", vertex_t::tag_t::infinity);
    EXPECT_EQ(2, v1.index);
    EXPECT_NE(std::string(), v1.data);
    EXPECT_EQ(std::string("1h"), v1.data);
    EXPECT_EQ(vertex_t::tag_t::infinity, v1.tag);
    EXPECT_TRUE(v1.is_infinity());

    vertex_t v2(4, "", vertex_t::tag_t::repeat);
    EXPECT_EQ(4, v2.index);
    EXPECT_EQ(std::string(), v2.data);
    EXPECT_EQ(vertex_t::tag_t::repeat, v2.tag);
    EXPECT_TRUE(!v2.is_infinity());
}

TEST(VertexTest, CompareFunction) {
    std::set<vertex_ptr, VertexPtrComporator<vertex_t>> verteces;

    verteces.insert(vertex_ptr(new vertex_t(1, "1h", vertex_t::tag_t::block)));
    verteces.insert(vertex_ptr(new vertex_t(2, "1t", vertex_t::tag_t::infinity)));
    verteces.insert(vertex_ptr(new vertex_t(4, "repeat", vertex_t::tag_t::repeat)));

    EXPECT_EQ((*verteces.begin())->index, 1);
    EXPECT_EQ((*(++verteces.begin()))->index, 2);
    EXPECT_EQ((*verteces.rbegin())->index, 4);
    EXPECT_EQ(verteces.size(), 3);

    verteces.insert(vertex_ptr(new vertex_t(1, "2h", vertex_t::tag_t::block)));
    EXPECT_EQ((*verteces.begin())->index, 1);
    EXPECT_EQ((*verteces.begin())->data, "1h");
    EXPECT_EQ(verteces.size(), 3);

    verteces.insert(vertex_ptr(new vertex_t(3, "2h", vertex_t::tag_t::block)));
    EXPECT_EQ((*(++verteces.rbegin()))->index, 3);
    EXPECT_EQ((*verteces.rbegin())->index, 4);
    EXPECT_EQ(verteces.size(), 4);

}