//
// Created by pavel on 11/6/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"
#include "graph/vertex.hpp"

using namespace std;

using vertex_t = Vertex<std::string>;

TEST(VertexTest, ConstructorWithOneParameter) {
    vertex_t v(1);
    EXPECT_EQ(1, v.index);
    EXPECT_EQ(std::string(), v.data);
    EXPECT_EQ(vertex_t::tag_t::block, v.tag);
    EXPECT_TRUE(!v.is_infinity());
    EXPECT_EQ(v, v);
}

TEST(VertexTest, ConstructorWithOneParamterPtr) {
    using vertex_ptr = vertex_t::vertex_ptr;

    vertex_ptr v_ptr(new vertex_t(1, "1t"));
    vertex_t v(v_ptr);
    EXPECT_EQ(1, v.index);
    EXPECT_EQ(std::string("1t"), v.data);
    EXPECT_EQ(vertex_t::tag_t::block, v.tag);
    EXPECT_TRUE(!v.is_infinity());
    EXPECT_EQ(v, *v_ptr);
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

    vertex_t v1(0, "1h", vertex_t::tag_t::infinity);
    EXPECT_EQ(0, v1.index);
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