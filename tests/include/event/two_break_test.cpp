//
// Created by pavel on 11/6/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace std;
using namespace event;
using namespace structure;
using vertex_t = Vertex<std::string>;
using mcolor_t = Mcolor;
using vertex_ptr = typename vertex_t::vertex_ptr;

class TwoBreakTest : public testing::Test {
    virtual void SetUp() override {
        verteces.push_back(vertex_ptr(new vertex_t(1, "1h")));
        verteces.push_back(vertex_ptr(new vertex_t(2, "1t")));
        verteces.push_back(vertex_ptr(new vertex_t(3, "2h")));
        verteces.push_back(vertex_ptr(new vertex_t(4, "2t")));
        verteces.push_back(vertex_ptr(new vertex_t(5, "3t")));
        verteces.push_back(vertex_ptr(new vertex_t(6, "3h")));
        verteces.push_back(vertex_ptr(new vertex_t(7, "oo", vertex_t::tag_t::infinity)));

        color_twobreaks.push_back(color_twobreak_t(verteces[0], verteces[1], verteces[2], verteces[3], Mcolor(1)));

        twobreaks.push_back(twobreak_t(verteces[0], verteces[1], verteces[2], verteces[3])); //0
        twobreaks.push_back(twobreak_t(std::make_pair(verteces[0], verteces[2]), std::make_pair(verteces[1], verteces[3]))); //1
        twobreaks.push_back(twobreak_t(color_twobreaks[0])); //2
        twobreaks.push_back(twobreak_t(verteces[0], verteces[1], verteces[3], verteces[2])); //3
        twobreaks.push_back(twobreak_t(verteces[1], verteces[0], verteces[2], verteces[3])); //4
        twobreaks.push_back(twobreak_t(verteces[1], verteces[0], verteces[3], verteces[2])); //5
        twobreaks.push_back(twobreak_t(verteces[2], verteces[3], verteces[0], verteces[1])); //6
        twobreaks.push_back(twobreak_t(verteces[2], verteces[3], verteces[1], verteces[0])); //7
        twobreaks.push_back(twobreak_t(verteces[3], verteces[2], verteces[0], verteces[1])); //8
        twobreaks.push_back(twobreak_t(verteces[3], verteces[2], verteces[1], verteces[0])); //9

        twobreaks.push_back(twobreak_t(verteces[0], verteces[1], verteces[5], verteces[6])); //10
        twobreaks.push_back(twobreak_t(verteces[1], verteces[0], verteces[5], verteces[6])); //11
        twobreaks.push_back(twobreak_t(verteces[1], verteces[0], verteces[6], verteces[5])); //12
        twobreaks.push_back(twobreak_t(verteces[5], verteces[6], verteces[0], verteces[1])); //13
        twobreaks.push_back(twobreak_t(verteces[5], verteces[6], verteces[1], verteces[0])); //14
        twobreaks.push_back(twobreak_t(verteces[6], verteces[5], verteces[0], verteces[1])); //15
        twobreaks.push_back(twobreak_t(verteces[6], verteces[5], verteces[1], verteces[0])); //16

        //TODO CONTINUE change
        //TODO ADD INFINITY
    }

protected:
    using twobreak_t = TwoBreak<vertex_t, void>;
    using color_twobreak_t = TwoBreak<vertex_t, Mcolor>;

    std::vector<vertex_ptr> verteces;
    std::vector<twobreak_t> twobreaks;
    std::vector<twobreak_t> color_twobreaks;
    //color_twobreak_t color_break;
};

TEST_F(TwoBreakTest, VertexConstructor) {
    EXPECT_EQ(*verteces[0], *(twobreaks[0].get_vertex(0)));
    EXPECT_EQ(*verteces[1], *(twobreaks[0].get_vertex(1)));
    EXPECT_EQ(*verteces[2], *(twobreaks[0].get_vertex(2)));
    EXPECT_EQ(*verteces[3], *(twobreaks[0].get_vertex(3)));
    EXPECT_EQ(*verteces[0], *(twobreaks[0].get_edge(0).first));
    EXPECT_EQ(*verteces[1], *(twobreaks[0].get_edge(0).second));
    EXPECT_EQ(*verteces[2], *(twobreaks[0].get_edge(1).first));
    EXPECT_EQ(*verteces[3], *(twobreaks[0].get_edge(1).second));
}

TEST_F(TwoBreakTest, EdgeConstructor) {
    EXPECT_EQ(*verteces[0], *(twobreaks[1].get_vertex(0)));
    EXPECT_EQ(*verteces[2], *(twobreaks[1].get_vertex(1)));
    EXPECT_EQ(*verteces[1], *(twobreaks[1].get_vertex(2)));
    EXPECT_EQ(*verteces[3], *(twobreaks[1].get_vertex(3)));
    EXPECT_EQ(*verteces[0], *(twobreaks[1].get_edge(0).first));
    EXPECT_EQ(*verteces[2], *(twobreaks[1].get_edge(0).second));
    EXPECT_EQ(*verteces[1], *(twobreaks[1].get_edge(1).first));
    EXPECT_EQ(*verteces[3], *(twobreaks[1].get_edge(1).second));
}

TEST_F(TwoBreakTest, ColorTwoBreakConstructor) {
    EXPECT_EQ(*verteces[0], *(twobreaks[2].get_vertex(0)));
    EXPECT_EQ(*verteces[1], *(twobreaks[2].get_vertex(1)));
    EXPECT_EQ(*verteces[2], *(twobreaks[2].get_vertex(2)));
    EXPECT_EQ(*verteces[3], *(twobreaks[2].get_vertex(3)));
    EXPECT_EQ(*verteces[0], *(twobreaks[2].get_edge(0).first));
    EXPECT_EQ(*verteces[1], *(twobreaks[2].get_edge(0).second));
    EXPECT_EQ(*verteces[2], *(twobreaks[2].get_edge(1).first));
    EXPECT_EQ(*verteces[3], *(twobreaks[2].get_edge(1).second));
}

TEST_F(TwoBreakTest, ChangeVertexFunction) {
    twobreaks[0].change_vertex(0, verteces[4]);
    twobreaks[0].change_vertex(1, verteces[5]);
    twobreaks[0].change_vertex(2, verteces[0]);
    twobreaks[0].change_vertex(3, verteces[1]);

    EXPECT_EQ(*verteces[4], *(twobreaks[0].get_vertex(0)));
    EXPECT_EQ(*verteces[5], *(twobreaks[0].get_vertex(1)));
    EXPECT_EQ(*verteces[0], *(twobreaks[0].get_vertex(2)));
    EXPECT_EQ(*verteces[1], *(twobreaks[0].get_vertex(3)));
}

TEST_F(TwoBreakTest, TestInverseFunction) {
    EXPECT_EQ(twobreaks[0], twobreaks[1].inverse());
    EXPECT_EQ(twobreaks[1], twobreaks[2].inverse());
}

TEST_F(TwoBreakTest, TestDependentFunction) {
    EXPECT_EQ(twobreaks[0].is_dependent(twobreaks[1]), twobreak_t::dependence_type::strong_dependent);
    EXPECT_EQ(twobreaks[5].is_dependent(twobreaks[1]), twobreak_t::dependence_type::strong_dependent);
    EXPECT_EQ(twobreaks[6].is_dependent(twobreaks[1]), twobreak_t::dependence_type::strong_dependent);
    EXPECT_EQ(twobreaks[9].is_dependent(twobreaks[1]), twobreak_t::dependence_type::strong_dependent);

    EXPECT_EQ(twobreaks[3].is_dependent(twobreaks[1]), twobreak_t::dependence_type::independent);
    EXPECT_EQ(twobreaks[4].is_dependent(twobreaks[1]), twobreak_t::dependence_type::independent);
    EXPECT_EQ(twobreaks[7].is_dependent(twobreaks[1]), twobreak_t::dependence_type::independent);
    EXPECT_EQ(twobreaks[8].is_dependent(twobreaks[1]), twobreak_t::dependence_type::independent);

    EXPECT_EQ(twobreaks[1].is_dependent(twobreaks[10]), twobreak_t::dependence_type::weakly_dependent);
    EXPECT_EQ(twobreaks[1].is_dependent(twobreaks[12]), twobreak_t::dependence_type::weakly_dependent);
    EXPECT_EQ(twobreaks[1].is_dependent(twobreaks[13]), twobreak_t::dependence_type::weakly_dependent);
    EXPECT_EQ(twobreaks[1].is_dependent(twobreaks[16]), twobreak_t::dependence_type::weakly_dependent);

    //TODO think many tests
}

/*
TEST(MulticolorTwoBreakTest, BasicConstructor) {
    using twobreak_t = TwoBreak<vertex_t, mcolor_t >;

    vertex_ptr v1(new vertex_t(1, "1h"));
    vertex_ptr v2(new vertex_t(2, "1t"));
    vertex_ptr v3(new vertex_t(3, "2h"));
    vertex_ptr v4(new vertex_t(4, "2t"));
    mcolor_t color(1);

    twobreak_t tbr(v1, v2, v3, v4, color);
    EXPECT_EQ(v1, tbr.get_vertex(0));
    EXPECT_EQ(v2, tbr.get_vertex(1));
    EXPECT_EQ(v3, tbr.get_vertex(2));
    EXPECT_EQ(v4, tbr.get_vertex(3));
    EXPECT_EQ(vertex_t(1, "1h"), *(tbr.get_vertex(0)));
    EXPECT_EQ(vertex_t(2, "1t"), *(tbr.get_vertex(1)));
    EXPECT_EQ(vertex_t(3, "2h"), *(tbr.get_vertex(2)));
    EXPECT_EQ(vertex_t(4, "2t"), *(tbr.get_vertex(3)));
    EXPECT_EQ(color, tbr.get_multicolor());
}*/
