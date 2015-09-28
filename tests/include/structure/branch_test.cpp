//
// Created by Nikita Kartashov on 30/03/2015.
//

#include "gtest/gtest.h"
#include "structures/branch.hpp"
#include "structures/mcolor.hpp"

using namespace structure;
using namespace std;

TEST(BranchTest, DoIntersect) {
  using BranchHelper = Branch<Mcolor>;
  using branch_t = pair<Mcolor, Mcolor>;
  Mcolor left(1);
  Mcolor right(2);
  auto left_branch = branch_t(left, right);
  auto right_branch = branch_t(Mcolor(left, right, Mcolor::Union), Mcolor());
  ASSERT_FALSE(BranchHelper::do_intersect(left_branch, right_branch));
}