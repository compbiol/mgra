//
// Created by pavel on 10/26/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace std;
using namespace event;
using namespace structure;
using wgd_t = wgd<std::string>;

TEST(WGDTest, BasicWgd) {
    wgd_t wgd("A", "B", std::vector<std::string>({"WGD1"}));
    EXPECT_EQ("A", wgd.get_parent());
    EXPECT_EQ("B", wgd.get_children());
    EXPECT_EQ(*wgd.cbegin(), *wgd.crbegin());
    EXPECT_EQ(*wgd.begin(), *wgd.crbegin());
    EXPECT_EQ(*wgd.cbegin(), "WGD1");
}

TEST(WGDTest, MultiplyWgd) {
    wgd_t wgd("A", "B", std::vector<std::string>({"WGD1", "WGD2", "WGD3"}));
    EXPECT_EQ("A", wgd.get_parent());
    EXPECT_EQ("B", wgd.get_children());
    EXPECT_EQ(*wgd.cbegin(), "WGD1");
    EXPECT_EQ(*wgd.begin(), "WGD1");
    EXPECT_EQ(*(--wgd.crend()), "WGD1");

    EXPECT_EQ(*(++wgd.cbegin()), *(++wgd.crbegin()));
    EXPECT_EQ(*(++wgd.cbegin()), "WGD2");
    EXPECT_EQ(*(++wgd.crbegin()), "WGD2");
    EXPECT_EQ(*(++wgd.begin()), "WGD2");

    EXPECT_EQ(*(--wgd.cend()), "WGD3");
    EXPECT_EQ(*(--wgd.end()), "WGD3");
    EXPECT_EQ(*wgd.crbegin(), "WGD3");
}