//
// Created by pavel on 12/28/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace structure;
using namespace structure::phyl_tree;

using string_tree_t = BinaryTree<std::string>;
using string_tree_builder_t = TreeBuilder<string_tree_t>;
using string_node_t = string_tree_t::node_t;
using string_node_ptr = string_tree_t::node_ptr;

using mcolor_tree_t = BinaryTree<Mcolor>;
using mcolor_tree_builder_t = TreeBuilder<mcolor_tree_t>;
using mcolor_node_t = mcolor_tree_t::node_t;
using mcolor_node_ptr = mcolor_tree_t::node_ptr;

//TODO here think about test more clearly

TEST(TreeBuilderTest, NewickBuildForOneGenomes) {
    std::vector<std::string> strs({"A", " A ", " A", "Rat ", " Rat ", " rat_1 "});

    for (auto const &str : strs) {
        string_tree_t tree = string_tree_builder_t::build_tree(str);

        EXPECT_EQ(tree.get_root()->get_parent(), nullptr);
        EXPECT_EQ(tree.get_root()->get_type(), string_node_t::classic);
        EXPECT_FALSE(tree.get_root()->has_childrens());
        EXPECT_TRUE(tree.get_root()->is_leaf());
    }

    for (size_t i = 0; i < 3; ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(strs[i]);
        EXPECT_EQ(tree.get_root()->get_data(), "A");
    }

    for (size_t i = 3; i < 5; ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(strs[i]);
        EXPECT_EQ(tree.get_root()->get_data(), "Rat");
    }

    string_tree_t tree = string_tree_builder_t::build_tree(strs[5]);
    EXPECT_EQ(tree.get_root()->get_data(), "rat_1");
}


TEST(TreeBuilderTest, NewickBuildForTwoGenomes) {
    std::vector<std::string> strs(
            {"(A,B)", "( A , B )", "( A,B)", "(A ,B)", "(A, B)", "(A,B )", "(A ,B):3", "(A, B) : 3 ", "(Dog,Rat)", "( Dog, Rat )", "(Dog, Rat)",
             "(Dog , Rat) : 3 ", "(Dog, Rat)DR", "(Dog, Rat) DR ", "(Dog , Rat) DR : 3 ", " ( Dog , Rat ) DR : 3 ", "(Dog,Rat)DR:3"});

    for (auto const & str : strs) {
        string_tree_t tree = string_tree_builder_t::build_tree(str);

        EXPECT_EQ(tree.get_root()->get_parent(), nullptr);
        EXPECT_EQ(tree.get_root()->get_type(), string_node_t::classic);
        EXPECT_TRUE(tree.get_root()->has_childrens());
        EXPECT_FALSE(tree.get_root()->is_leaf());

        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_parent(), tree.get_root().get());
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_type(), string_node_t::classic);
        EXPECT_TRUE(tree.get_root()->get_most_left_child()->is_leaf());
        EXPECT_FALSE(tree.get_root()->get_most_left_child()->has_childrens());

        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_parent(), tree.get_root().get());
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_type(), string_node_t::classic);
        EXPECT_TRUE(tree.get_root()->get_most_right_child()->is_leaf());
        EXPECT_FALSE(tree.get_root()->get_most_right_child()->has_childrens());
    }

    for (size_t i = 0; i < 8; ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(strs[i]);
        EXPECT_EQ(tree.get_root()->get_data(), "");
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "A");
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_data(), "B");
    }

    for (size_t i = 8; i < 12; ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(strs[i]);
        EXPECT_EQ(tree.get_root()->get_data(), "");
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "Dog");
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_data(), "Rat");
    }

    for (size_t i = 12; i < strs.size(); ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(strs[i]);
        EXPECT_EQ(tree.get_root()->get_data(), "DR");
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "Dog");
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_data(), "Rat");
    }
}

TEST(TreeBuilderTest, NewickBuildBinaryTreeForThreeGenomes) {
    std::vector<std::string> first_strs({"((A,B)AB:3,C) ", "( ( A , B ) AB : 3, C ) ", "((A,B)AB,C)"});

    for (auto const &str : first_strs) {
        string_tree_t tree = string_tree_builder_t::build_tree(str);

        EXPECT_EQ(tree.get_root()->get_parent(), nullptr);
        EXPECT_EQ(tree.get_root()->get_data(), "");
        EXPECT_EQ(tree.get_root()->get_type(), string_node_t::classic);
        EXPECT_TRUE(tree.get_root()->has_childrens());
        EXPECT_FALSE(tree.get_root()->is_leaf());

        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "AB");
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_parent(), tree.get_root().get());
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_type(), string_node_t::classic);
        EXPECT_TRUE(tree.get_root()->get_most_left_child()->has_childrens());
        EXPECT_FALSE(tree.get_root()->get_most_left_child()->is_leaf());

        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_data(), "C");
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_parent(), tree.get_root().get());
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_type(), string_node_t::classic);
        EXPECT_FALSE(tree.get_root()->get_most_right_child()->has_childrens());
        EXPECT_TRUE(tree.get_root()->get_most_right_child()->is_leaf());
    }

    std::vector<std::string> second_strs(
            {"(Cat,(Rat,Dog)RD:3)", "(  Cat, ( Rat , Dog ) RD : 3 )", "(Cat,(Rat,Dog)RD)"});
    for (auto const &str : second_strs) {
        string_tree_t tree = string_tree_builder_t::build_tree(str);

        EXPECT_EQ(tree.get_root()->get_parent(), nullptr);
        EXPECT_EQ(tree.get_root()->get_data(), "");
        EXPECT_EQ(tree.get_root()->get_type(), string_node_t::classic);
        EXPECT_TRUE(tree.get_root()->has_childrens());
        EXPECT_FALSE(tree.get_root()->is_leaf());

        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "Cat");
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_parent(), tree.get_root().get());
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_type(), string_node_t::classic);
        EXPECT_FALSE(tree.get_root()->get_most_left_child()->has_childrens());
        EXPECT_TRUE(tree.get_root()->get_most_left_child()->is_leaf());

        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_data(), "RD");
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_parent(), tree.get_root().get());
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_type(), string_node_t::classic);
        EXPECT_TRUE(tree.get_root()->get_most_right_child()->has_childrens());
        EXPECT_FALSE(tree.get_root()->get_most_right_child()->is_leaf());
    }
}

TEST(TreeBuilderTest, NewickBuildGeneralTreeForThreeGenomes) {
    std::vector<std::string> third_strs({"(A,B,C)ABC", " ( A , B , C ) ABC ", "(A, B, C) ABC : 4"});
    for (auto const &str : third_strs) {
        string_tree_t tree = string_tree_builder_t::build_tree(str);

        EXPECT_EQ(tree.get_root()->get_parent(), nullptr);
        EXPECT_EQ(tree.get_root()->get_data(), "ABC");
        EXPECT_EQ(tree.get_root()->get_type(), string_node_t::subtree);
        EXPECT_TRUE(tree.get_root()->has_childrens());
        EXPECT_FALSE(tree.get_root()->is_leaf());

        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "A");
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_parent(), tree.get_root().get());
        EXPECT_EQ(tree.get_root()->get_most_left_child()->get_type(), string_node_t::classic);
        EXPECT_FALSE(tree.get_root()->get_most_left_child()->has_childrens());
        EXPECT_TRUE(tree.get_root()->get_most_left_child()->is_leaf());

        EXPECT_EQ((*++tree.get_root()->cbegin())->get_data(), "B");
        EXPECT_EQ((*++tree.get_root()->cbegin())->get_parent(), tree.get_root().get());
        EXPECT_EQ((*++tree.get_root()->cbegin())->get_type(), string_node_t::classic);
        EXPECT_FALSE((*++tree.get_root()->cbegin())->has_childrens());
        EXPECT_TRUE((*++tree.get_root()->cbegin())->is_leaf());

        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_data(), "C");
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_parent(), tree.get_root().get());
        EXPECT_EQ(tree.get_root()->get_most_right_child()->get_type(), string_node_t::classic);
        EXPECT_FALSE(tree.get_root()->get_most_right_child()->has_childrens());
        EXPECT_TRUE(tree.get_root()->get_most_right_child()->is_leaf());
    }
}

TEST(TreeBuilderTest, NewickBuildBinaryForNGenomes) {
    //TODO here write binary tree
}

TEST(TreeBuilderTest, NewickBuildForNGenomes) {
    std::string str = "((((A,B,C)ABC:3,(D,H,Mouse)DHM:4)ABCDHM:5,(Rat,Dog)FE)ABCDHMFE:5,((Cat, (Human, Opossum)3)4,Chimpanzee)5:3)";

    string_tree_t tree = string_tree_builder_t::build_tree(str);

    EXPECT_EQ(tree.get_root()->get_parent(), nullptr);
    EXPECT_EQ(tree.get_root()->get_data(), "");
    EXPECT_EQ(tree.get_root()->get_type(), string_node_t::classic);
    EXPECT_TRUE(tree.get_root()->has_childrens());
    EXPECT_FALSE(tree.get_root()->is_leaf());

    auto first_left_child = tree.get_root()->get_most_left_child();
    EXPECT_EQ(first_left_child->get_data(), "ABCDHMFE");
    EXPECT_EQ(first_left_child->get_parent(), tree.get_root().get());
    EXPECT_EQ(first_left_child->get_type(), string_node_t::classic);
    EXPECT_TRUE(first_left_child->has_childrens());
    EXPECT_FALSE(first_left_child->is_leaf());

    auto first_right_child = tree.get_root()->get_most_right_child();
    EXPECT_EQ(first_right_child->get_data(), "5");
    EXPECT_EQ(first_right_child->get_parent(), tree.get_root().get());
    EXPECT_EQ(first_right_child->get_type(), string_node_t::classic);
    EXPECT_TRUE(first_right_child->has_childrens());
    EXPECT_FALSE(first_right_child->is_leaf());

    //TODO here add checks
}


//TODO TODO branches and str for n genomes