//
// Created by Nikita Kartashov on 30/03/2015.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace structure;
using branch_t = Branch<Mcolor>;

TEST(BranchTest, BasicConstructor) {
    Mcolor left(1); left.insert(2); left.insert(3);
    Mcolor right(4);
    branch_t branch1(left, right);
    branch_t branch2(right, left);
    EXPECT_EQ(branch1.left, left);
    EXPECT_EQ(branch1.right, right);
    EXPECT_EQ(branch2.left, left);
    EXPECT_EQ(branch2.right, right);
}

TEST(BranchTest, DoIntersect) {
    Mcolor left(1); Mcolor right(2);
    auto left_branch = branch_t(left, right);
    auto right_branch = branch_t(Mcolor(left, right, Mcolor::Union), Mcolor());
    ASSERT_FALSE(left_branch.do_intersect(right_branch));

    //TODO think many tests
}