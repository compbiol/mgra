//
// Created by pavel on 12/28/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace structure;
using namespace structure::phyl_tree;

using color_tree_t = BinaryTree<Mcolor>;
using color_node_t = color_tree_t::node_t;
using color_node_ptr = color_tree_t::node_ptr;

using string_tree_t = BinaryTree<std::string>;
using string_node_t = string_tree_t::node_t;
using string_node_ptr = string_tree_t::node_ptr;

TEST(TreeTest, DefaultConstructiorColorTree) {
    color_tree_t tree;
    EXPECT_EQ(tree.get_root(), nullptr);
}

TEST(TreeTest, BasicConstructiorColorTree) {
    color_node_ptr root(new color_node_t(Mcolor({1,2}), color_node_ptr(new color_node_t(Mcolor(1))), color_node_ptr(new color_node_t(Mcolor(2)))));
    color_tree_t tree(root);
    EXPECT_EQ(tree.get_root(), root);
}

TEST(TreeTest, DefaultConstructiorStringTree) {
    string_tree_t tree;
    EXPECT_EQ(tree.get_root(), nullptr);
}

TEST(TreeTest, BasicConstructiorStringTree) {
    string_node_ptr root(new string_node_t("AB", string_node_ptr(new string_node_t("A")), string_node_ptr(new string_node_t("B"))));
    string_tree_t tree(root);
    EXPECT_EQ(tree.get_root(), root);
}

