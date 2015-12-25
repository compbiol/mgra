//
// Created by pavel on 12/23/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"
#include "structures/mcolor.hpp"

using namespace structure;
using namespace structure::phyl_tree;

using node_t = Node<Mcolor>;
using node_ptr = std::shared_ptr<Node<Mcolor>>;

TEST(NodeTest, ClassicLeaf) {
    node_t node(Mcolor(1));
    EXPECT_EQ(node.get_data(), Mcolor(1));
    EXPECT_EQ(node.get_type(), node_t::classic);
    EXPECT_EQ(node.get_left_child(), nullptr);
    EXPECT_EQ(node.get_right_child(), nullptr);
    EXPECT_EQ(node.get_parent(), nullptr);
    EXPECT_FALSE(node.has_left_child());
    EXPECT_FALSE(node.has_right_child());
    EXPECT_TRUE(node.is_leaf());
    EXPECT_TRUE(node.is_basic_leaf());
    EXPECT_FALSE(node.is_subtree_leaf());
}

TEST(NodeTest, SubtreeLeaf) {
    node_t node(Mcolor({1, 2, 3}));
    EXPECT_EQ(node.get_data(), Mcolor({1, 2, 3}));
    EXPECT_EQ(node.get_type(), node_t::classic);
    EXPECT_EQ(node.get_left_child(), nullptr);
    EXPECT_EQ(node.get_right_child(), nullptr);
    EXPECT_EQ(node.get_parent(), nullptr);
    EXPECT_FALSE(node.has_left_child());
    EXPECT_FALSE(node.has_right_child());
    EXPECT_TRUE(node.is_leaf());
    EXPECT_FALSE(node.is_basic_leaf());
    EXPECT_TRUE(node.is_subtree_leaf());
}

TEST(NodeTest, InternalClassicNode) {
    node_ptr child_1(new node_t(Mcolor(1)));
    node_ptr child_2(new node_t(Mcolor({2, 3})));
    node_ptr internal_node(new node_t(Mcolor({1, 2, 3}), child_1, child_2));
    child_1->set_parent(internal_node.get());
    child_2->set_parent(internal_node.get());

    EXPECT_EQ(internal_node->get_data(), Mcolor({1, 2, 3}));
    EXPECT_EQ(internal_node->get_left_child(), child_1);
    EXPECT_EQ(internal_node->get_right_child(), child_2);
    EXPECT_EQ(internal_node->get_type(), node_t::classic);
    EXPECT_EQ(internal_node->get_parent(), nullptr);
    EXPECT_TRUE(internal_node->has_left_child());
    EXPECT_TRUE(internal_node->has_right_child());
    EXPECT_FALSE(internal_node->is_leaf());
    EXPECT_FALSE(internal_node->is_basic_leaf());
    EXPECT_FALSE(internal_node->is_subtree_leaf());

    EXPECT_EQ(internal_node->get_left_child()->get_data(), Mcolor(1));
    EXPECT_EQ(internal_node->get_left_child()->get_parent(), internal_node.get());
    EXPECT_EQ(internal_node->get_left_child()->get_type(), node_t::classic);
    EXPECT_TRUE(internal_node->get_left_child()->is_leaf());
    EXPECT_TRUE(internal_node->get_left_child()->is_basic_leaf());
    EXPECT_FALSE(internal_node->get_left_child()->is_subtree_leaf());

    EXPECT_EQ(internal_node->get_right_child()->get_data(), Mcolor({2, 3}));
    EXPECT_EQ(internal_node->get_right_child()->get_parent(), internal_node.get());
    EXPECT_EQ(internal_node->get_right_child()->get_type(), node_t::classic);
    EXPECT_TRUE(internal_node->get_right_child()->is_leaf());
    EXPECT_FALSE(internal_node->get_right_child()->is_basic_leaf());
    EXPECT_TRUE(internal_node->get_right_child()->is_subtree_leaf());

    node_ptr leaf_1(new node_t(Mcolor(0)));
    node_ptr leaf_2(new node_t(Mcolor(2)));
    node_ptr leaf_3(new node_t(Mcolor(3)));

    internal_node->get_left_child()->set_left_child(leaf_1);
    internal_node->get_right_child()->set_left_child(leaf_2);
    internal_node->get_right_child()->set_right_child(leaf_2);

    EXPECT_FALSE(internal_node->get_left_child()->is_leaf());
    EXPECT_FALSE(internal_node->get_left_child()->is_basic_leaf());
    EXPECT_FALSE(internal_node->get_left_child()->is_subtree_leaf());

    EXPECT_FALSE(internal_node->get_right_child()->is_leaf());
    EXPECT_FALSE(internal_node->get_right_child()->is_basic_leaf());
    EXPECT_FALSE(internal_node->get_right_child()->is_subtree_leaf());

    //TODO Finish here, check get_parents
}


TEST(NodeTest, WgdNode) {
    node_ptr child(new node_t(Mcolor(1)));
    node_ptr wgd_node(new node_t(Mcolor({1, 1}), child));
    child->set_parent(wgd_node.get());

    EXPECT_EQ(wgd_node->get_data(), Mcolor({1, 1}));
    EXPECT_EQ(wgd_node->get_left_child(), child);
    EXPECT_EQ(wgd_node->get_right_child(), nullptr);
    EXPECT_EQ(wgd_node->get_type(), node_t::whole_duplication);
    EXPECT_EQ(wgd_node->get_parent(), nullptr);
    EXPECT_TRUE(wgd_node->has_left_child());
    EXPECT_FALSE(wgd_node->has_right_child());
    EXPECT_FALSE(wgd_node->is_leaf());
    EXPECT_FALSE(wgd_node->is_basic_leaf());
    EXPECT_FALSE(wgd_node->is_subtree_leaf());

    EXPECT_EQ(wgd_node->get_left_child()->get_data(), Mcolor(1));
    EXPECT_EQ(wgd_node->get_left_child()->get_parent(), wgd_node.get());
    EXPECT_EQ(wgd_node->get_left_child()->get_type(), node_t::classic);
}