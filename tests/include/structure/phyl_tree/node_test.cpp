//
// Created by pavel on 12/23/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace structure;
using namespace structure::phyl_tree;

using string_node_t = Node<std::string>;
using string_node_ptr = std::shared_ptr<string_node_t>;

TEST(NodeTest, ClassicLeaf) {
    string_node_t node("A");
    EXPECT_EQ(node.get_data(), "A");
    EXPECT_EQ(node.get_type(), string_node_t::classic);
    EXPECT_FALSE(node.has_childrens());
    EXPECT_EQ(node.get_parent(), nullptr);
    EXPECT_TRUE(node.is_leaf());
}

TEST(NodeTest, WgdNode) {
    string_node_ptr child(new string_node_t("A"));
    string_node_ptr wgd_node(new string_node_t("Wgd1", child));

    EXPECT_EQ(wgd_node->get_data(), "Wgd1");
    EXPECT_EQ(wgd_node->get_type(), string_node_t::whole_duplication);
    EXPECT_TRUE(wgd_node->has_childrens());
    EXPECT_EQ(wgd_node->get_most_left_child(), child);
    EXPECT_EQ(wgd_node->get_most_right_child(), child);
    EXPECT_EQ(wgd_node->get_parent(), nullptr);
    EXPECT_FALSE(wgd_node->is_leaf());

    EXPECT_EQ(wgd_node->get_most_left_child()->get_data(), "A");
    EXPECT_EQ(wgd_node->get_most_left_child()->get_type(), string_node_t::classic);
    EXPECT_EQ(wgd_node->get_most_left_child()->get_parent(), wgd_node.get());
    EXPECT_FALSE(wgd_node->get_most_left_child()->has_childrens());
    EXPECT_TRUE(wgd_node->get_most_left_child()->is_leaf());
}

TEST(NodeTest, ClassicNode) {
    string_node_ptr child1(new string_node_t("A"));
    string_node_ptr child2(new string_node_t("B"));
    string_node_ptr node(new string_node_t("AB", child1, child2));

    EXPECT_EQ(node->get_data(), "AB");
    EXPECT_EQ(node->get_type(), string_node_t::classic);
    EXPECT_TRUE(node->has_childrens());
    EXPECT_EQ(node->get_most_left_child(), child1);
    EXPECT_EQ(node->get_most_right_child(), child2);
    EXPECT_EQ(node->get_parent(), nullptr);
    EXPECT_FALSE(node->is_leaf());

    EXPECT_EQ(node->get_most_left_child()->get_data(), "A");
    EXPECT_EQ(node->get_most_left_child()->get_type(), string_node_t::classic);
    EXPECT_EQ(node->get_most_left_child()->get_parent(), node.get());
    EXPECT_FALSE(node->get_most_left_child()->has_childrens());
    EXPECT_TRUE(node->get_most_left_child()->is_leaf());

    EXPECT_EQ(node->get_most_right_child()->get_data(), "B");
    EXPECT_EQ(node->get_most_right_child()->get_type(), string_node_t::classic);
    EXPECT_EQ(node->get_most_right_child()->get_parent(), node.get());
    EXPECT_FALSE(node->get_most_right_child()->has_childrens());
    EXPECT_TRUE(node->get_most_right_child()->is_leaf());
}

TEST(NodeTest, SubtreeNode) {
    std::vector<string_node_ptr> childs;
    childs.push_back(string_node_ptr(new string_node_t("A")));
    childs.push_back(string_node_ptr(new string_node_t("B")));
    childs.push_back(string_node_ptr(new string_node_t("C")));

    string_node_ptr node(new string_node_t("ABC", childs));
    EXPECT_EQ(node->get_data(), "ABC");
    EXPECT_EQ(node->get_type(), string_node_t::subtree);
    EXPECT_TRUE(node->has_childrens());
    EXPECT_EQ(node->get_most_left_child(), childs[0]);
    EXPECT_EQ(*(++node->cbegin()), childs[1]);
    EXPECT_EQ(node->get_most_right_child(), childs[2]);
    EXPECT_EQ(node->get_parent(), nullptr);
    EXPECT_FALSE(node->is_leaf());

    EXPECT_EQ(node->get_most_left_child()->get_data(), "A");
    EXPECT_EQ(node->get_most_left_child()->get_type(), string_node_t::classic);
    EXPECT_EQ(node->get_most_left_child()->get_parent(), node.get());
    EXPECT_FALSE(node->get_most_left_child()->has_childrens());
    EXPECT_TRUE(node->get_most_left_child()->is_leaf());

    EXPECT_EQ((*(++node->cbegin()))->get_data(), "B");
    EXPECT_EQ((*(++node->cbegin()))->get_type(), string_node_t::classic);
    EXPECT_EQ((*(++node->cbegin()))->get_parent(), node.get());
    EXPECT_FALSE((*(++node->cbegin()))->has_childrens());
    EXPECT_TRUE((*(++node->cbegin()))->is_leaf());

    EXPECT_EQ(node->get_most_right_child()->get_data(), "C");
    EXPECT_EQ(node->get_most_right_child()->get_type(), string_node_t::classic);
    EXPECT_EQ(node->get_most_right_child()->get_parent(), node.get());
    EXPECT_FALSE(node->get_most_right_child()->has_childrens());
    EXPECT_TRUE(node->get_most_right_child()->is_leaf());
}

TEST(NodeTest, AddChildrenFunction) {
    string_node_ptr child_1(new string_node_t("A"));
    string_node_ptr child_2(new string_node_t("BC"));
    string_node_ptr internal_node(new string_node_t("ABC", child_1, child_2));

    EXPECT_EQ(internal_node->get_data(), "ABC");
    EXPECT_TRUE(internal_node->has_childrens());
    EXPECT_EQ(internal_node->get_most_left_child(), child_1);
    EXPECT_EQ(internal_node->get_most_right_child(), child_2);
    EXPECT_EQ(internal_node->get_type(), string_node_t::classic);
    EXPECT_EQ(internal_node->get_parent(), nullptr);
    EXPECT_FALSE(internal_node->is_leaf());

    EXPECT_EQ(internal_node->get_most_left_child()->get_data(), "A");
    EXPECT_EQ(internal_node->get_most_left_child()->get_parent(), internal_node.get());
    EXPECT_EQ(internal_node->get_most_left_child()->get_type(), string_node_t::classic);
    EXPECT_TRUE(internal_node->get_most_left_child()->is_leaf());

    EXPECT_EQ(internal_node->get_most_right_child()->get_data(), "BC");
    EXPECT_EQ(internal_node->get_most_right_child()->get_parent(), internal_node.get());
    EXPECT_EQ(internal_node->get_most_right_child()->get_type(), string_node_t::classic);
    EXPECT_TRUE(internal_node->get_most_right_child()->is_leaf());

    string_node_ptr leaf_1(new string_node_t("B"));
    string_node_ptr leaf_2(new string_node_t("C"));
    string_node_ptr leaf_3(new string_node_t("D"));

    //TEST SET FUNCTIONS
    internal_node->get_most_right_child()->add_children(leaf_1);
    internal_node->get_most_right_child()->add_children(leaf_2);

    EXPECT_TRUE(internal_node->get_most_left_child()->is_leaf());
    EXPECT_FALSE(internal_node->get_most_left_child()->has_childrens());
    EXPECT_FALSE(internal_node->get_most_right_child()->is_leaf());
    EXPECT_TRUE(internal_node->get_most_right_child()->has_childrens());

    EXPECT_TRUE(internal_node->get_most_right_child()->get_most_left_child()->is_leaf());
    EXPECT_FALSE(internal_node->get_most_right_child()->get_most_left_child()->has_childrens());
    EXPECT_EQ(internal_node->get_most_right_child()->get_most_left_child()->get_parent(), internal_node->get_most_right_child().get());
    EXPECT_EQ(internal_node->get_most_right_child()->get_most_left_child()->get_type(), string_node_t::classic);

    EXPECT_TRUE(internal_node->get_most_right_child()->get_most_right_child()->is_leaf());
    EXPECT_FALSE(internal_node->get_most_right_child()->get_most_right_child()->has_childrens());
    EXPECT_EQ(internal_node->get_most_right_child()->get_most_right_child()->get_parent(), internal_node->get_most_right_child().get());
    EXPECT_EQ(internal_node->get_most_right_child()->get_most_right_child()->get_type(), string_node_t::classic);

    internal_node->add_children(leaf_3);
    EXPECT_EQ(internal_node->get_data(), "ABC");
    EXPECT_TRUE(internal_node->has_childrens());
    EXPECT_EQ(internal_node->get_type(), string_node_t::subtree);
    EXPECT_EQ(internal_node->get_parent(), nullptr);
    EXPECT_FALSE(internal_node->is_leaf());

    EXPECT_TRUE(internal_node->get_most_left_child()->is_leaf());
    EXPECT_FALSE(internal_node->get_most_left_child()->has_childrens());
    EXPECT_FALSE((*(++internal_node->cbegin()))->is_leaf());
    EXPECT_TRUE((*(++internal_node->cbegin()))->has_childrens());
    EXPECT_TRUE(internal_node->get_most_right_child()->is_leaf());
    EXPECT_FALSE(internal_node->get_most_right_child()->has_childrens());
}

TEST(NodeTest, ReplaceChildrenFunction) {
    string_node_ptr child_1(new string_node_t("A"));
    string_node_ptr child_2(new string_node_t("BC"));
    string_node_ptr internal_node(new string_node_t("ABC", child_1, child_2));

    EXPECT_EQ(internal_node->get_data(), "ABC");
    EXPECT_TRUE(internal_node->has_childrens());
    EXPECT_EQ(internal_node->get_most_left_child(), child_1);
    EXPECT_EQ(internal_node->get_most_right_child(), child_2);
    EXPECT_EQ(internal_node->get_type(), string_node_t::classic);
    EXPECT_EQ(internal_node->get_parent(), nullptr);
    EXPECT_FALSE(internal_node->is_leaf());

    EXPECT_EQ(internal_node->get_most_left_child()->get_data(), "A");
    EXPECT_EQ(internal_node->get_most_left_child()->get_parent(), internal_node.get());
    EXPECT_EQ(internal_node->get_most_left_child()->get_type(), string_node_t::classic);
    EXPECT_TRUE(internal_node->get_most_left_child()->is_leaf());

    EXPECT_EQ(internal_node->get_most_right_child()->get_data(), "BC");
    EXPECT_EQ(internal_node->get_most_right_child()->get_parent(), internal_node.get());
    EXPECT_EQ(internal_node->get_most_right_child()->get_type(), string_node_t::classic);
    EXPECT_TRUE(internal_node->get_most_right_child()->is_leaf());

    string_node_ptr child_3(new string_node_t("B"));
    string_node_ptr child_4(new string_node_t("C"));

    internal_node->replace_children(child_3, 0);
    internal_node->replace_children(child_4, 1);

    EXPECT_EQ(internal_node->get_data(), "ABC");
    EXPECT_TRUE(internal_node->has_childrens());
    EXPECT_EQ(internal_node->get_most_left_child(), child_3);
    EXPECT_EQ(internal_node->get_most_right_child(), child_4);
    EXPECT_EQ(internal_node->get_type(), string_node_t::classic);
    EXPECT_EQ(internal_node->get_parent(), nullptr);
    EXPECT_FALSE(internal_node->is_leaf());

    EXPECT_EQ(internal_node->get_most_left_child()->get_data(), "B");
    EXPECT_EQ(internal_node->get_most_left_child()->get_parent(), internal_node.get());
    EXPECT_EQ(internal_node->get_most_left_child()->get_type(), string_node_t::classic);
    EXPECT_TRUE(internal_node->get_most_left_child()->is_leaf());

    EXPECT_EQ(internal_node->get_most_right_child()->get_data(), "C");
    EXPECT_EQ(internal_node->get_most_right_child()->get_parent(), internal_node.get());
    EXPECT_EQ(internal_node->get_most_right_child()->get_type(), string_node_t::classic);
    EXPECT_TRUE(internal_node->get_most_right_child()->is_leaf());
}