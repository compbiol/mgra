//
// Created by Nikita Kartashov on 30/03/2015.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace structure;
using namespace structure::phyl_tree;
using node_t = Node<Mcolor>;
using branch_t = Branch<Mcolor>;

TEST(BranchTest, BasicConstructor) {
    Mcolor left({1, 2, 3});
    Mcolor right(4);
    branch_t branch1(left, right);
    branch_t branch2(right, left);
    EXPECT_EQ(branch1.left, left);
    EXPECT_EQ(branch1.right, right);
    EXPECT_EQ(branch2.left, left);
    EXPECT_EQ(branch2.right, right);
}

TEST(BranchTest, ConvertToNode) {
    Mcolor left({1, 1, 1, 2, 3});
    Mcolor right({4, 5, 6});

    branch_t branch(right, left);
    EXPECT_EQ(branch.left, left);
    EXPECT_EQ(branch.right, right);

    auto node = branch.convert_to_node();
    EXPECT_EQ(node.get_data(), Mcolor({1, 1, 1, 2, 3, 4, 5, 6}));
    EXPECT_EQ(node.get_parent(), nullptr);
    EXPECT_EQ(node.get_type(), node_t::classic);
    EXPECT_TRUE(node.has_childrens());
    EXPECT_FALSE(node.is_leaf());
    EXPECT_EQ(node.get_most_left_child()->get_data(), left);
    EXPECT_EQ(node.get_most_left_child()->get_type(), node_t::classic);
    EXPECT_TRUE(node.get_most_left_child()->is_leaf());
    EXPECT_FALSE(node.get_most_left_child()->has_childrens());
    EXPECT_EQ(node.get_most_right_child()->get_data(), right);
    EXPECT_EQ(node.get_most_right_child()->get_type(), node_t::classic);
    EXPECT_TRUE(node.get_most_right_child()->is_leaf());
    EXPECT_FALSE(node.get_most_right_child()->has_childrens());
}


TEST(BranchTest, DoIntersect) {
    /*Mcolor left(1); Mcolor right(2);
    auto left_branch = branch_t(left, right);
    auto right_branch = branch_t(Mcolor(left, right, Mcolor::Union), Mcolor());
    ASSERT_FALSE(left_branch.do_intersect(right_branch));*/

    //TODO think many tests
}