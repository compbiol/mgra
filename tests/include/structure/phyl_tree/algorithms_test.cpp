//
// Created by pavel on 12/29/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace event;
using namespace structure;
using namespace structure::phyl_tree;

using mcolor_t = Mcolor;
using wgd_t = event::wgd<std::string>;
using string_tree_t = BinaryTree<std::string>;
using node_t = string_tree_t::node_t;
using node_ptr = string_tree_t::node_ptr;
using string_tree_builder_t = TreeBuilder<string_tree_t>;

using mcolor_name_t = std::map<mcolor_t, std::string>;
using name_mcolor_t = std::unordered_map<std::string, mcolor_t>;
using medians_color_t = std::vector<std::tuple<mcolor_t, mcolor_t, mcolor_t>>;

//TODO organize test more clearly

class AlgorihmsTest : public testing::Test {
    virtual void SetUp() override {
        all_wgds.push_back(std::make_pair(std::make_pair("AB", "A"), wgd_t("AB", "A", std::vector<std::string>({"WGD1"}))));
        all_wgds.push_back(std::make_pair(std::make_pair("AB", "A"), wgd_t("AB", "A", std::vector<std::string>({"WGD1", "WGD2", "WGD3"}))));
        all_wgds.push_back(std::make_pair(std::make_pair("root", "AB"), wgd_t("root", "AB", std::vector<std::string>({"WGD2"}))));
        all_wgds.push_back(std::make_pair(std::make_pair("AB", "D"), wgd_t("AB", "D", std::vector<std::string>({"WGD4"}))));
        all_wgds.push_back(std::make_pair(std::make_pair("ABC", "C"), wgd_t("ABC", "C", std::vector<std::string>({"WGD5"}))));
        all_wgds.push_back(std::make_pair(std::make_pair("ABC", "AB"), wgd_t("ABC", "AB", std::vector<std::string>({"WGD6"}))));
        all_wgds.push_back(std::make_pair(std::make_pair("ABCD", "C"), wgd_t("ABCD", "C", std::vector<std::string>({"WGD7"}))));
        all_wgds.push_back(std::make_pair(std::make_pair("ABCD", "AB"), wgd_t("ABCD", "AB", std::vector<std::string>({"WGD8"}))));

        genome_number.insert(std::make_pair("A", 1));
        genome_number.insert(std::make_pair("B", 2));
        genome_number.insert(std::make_pair("C", 3));
        genome_number.insert(std::make_pair("D", 4));
        genome_number.insert(std::make_pair("E", 5));
        genome_number.insert(std::make_pair("F", 6));
    }

protected:
    using map_wgd_t = std::map<std::pair<std::string, std::string>, wgd_t>;
    using genome_number_t = std::unordered_map<std::string, size_t>;

    name_mcolor_t swipe(mcolor_name_t const & mcolor_name) const {
        name_mcolor_t result;
        for (auto const & info : mcolor_name) {
            assert(result.count(info.second) == 0);
            result.insert(std::make_pair(info.second, info.first));
        }
        return result;
    }

    std::vector<std::pair<std::pair<std::string, std::string>, wgd_t>> all_wgds;
    genome_number_t genome_number;
};

TEST_F(AlgorihmsTest, TwoGenomes) {
    map_wgd_t wgds({all_wgds[0], all_wgds[3]});

    string_tree_t tree = string_tree_builder_t::build_tree("(A,B)AB");
    refine_tree_by_wgds(tree, wgds);

    EXPECT_EQ(tree.get_root()->get_data(), "AB");
    EXPECT_EQ(tree.get_root()->get_type(), node_t::classic);
    EXPECT_TRUE(tree.get_root()->has_childrens());

    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "WGD1");
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_type(), node_t::whole_duplication);
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_parent(), tree.get_root().get());
    EXPECT_TRUE(tree.get_root()->get_most_left_child()->has_childrens());
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child(), tree.get_root()->get_most_left_child()->get_most_right_child());
    EXPECT_FALSE(tree.get_root()->get_most_left_child()->is_leaf());

    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child()->get_data(), "A");
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child()->get_type(), node_t::classic);
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child()->get_parent(), tree.get_root()->get_most_left_child().get());
    EXPECT_FALSE(tree.get_root()->get_most_left_child()->get_most_left_child()->has_childrens());
    EXPECT_TRUE(tree.get_root()->get_most_left_child()->get_most_left_child()->is_leaf());

    wgds = map_wgd_t({all_wgds[1]});
    tree = string_tree_builder_t::build_tree("(A,B)AB");
    refine_tree_by_wgds(tree, wgds);

    EXPECT_EQ(tree.get_root()->get_data(), "AB");
    EXPECT_EQ(tree.get_root()->get_type(), node_t::classic);
    EXPECT_TRUE(tree.get_root()->has_childrens());

    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "WGD1");
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_type(), node_t::whole_duplication);
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_parent(), tree.get_root().get());
    EXPECT_TRUE(tree.get_root()->get_most_left_child()->has_childrens());
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child(), tree.get_root()->get_most_left_child()->get_most_right_child());
    EXPECT_FALSE(tree.get_root()->get_most_left_child()->is_leaf());

    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child()->get_data(), "WGD2");
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child()->get_type(), node_t::whole_duplication);
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child()->get_parent(), tree.get_root()->get_most_left_child().get());
    EXPECT_TRUE(tree.get_root()->get_most_left_child()->get_most_left_child()->has_childrens());
    EXPECT_FALSE(tree.get_root()->get_most_left_child()->get_most_left_child()->is_leaf());
}

TEST_F(AlgorihmsTest, NGenomes){
    map_wgd_t wgds({all_wgds[1], all_wgds[2], all_wgds[3]});

    string_tree_t tree = string_tree_builder_t::build_tree("((A, B)AB, C, (D, F)DF, E)root");
    refine_tree_by_wgds(tree, wgds);

    EXPECT_EQ(tree.get_root()->get_data(), "root");
    EXPECT_EQ(tree.get_root()->get_type(), node_t::subtree);
    EXPECT_TRUE(tree.get_root()->has_childrens());

    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_data(), "WGD2");
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_type(), node_t::whole_duplication);
    EXPECT_TRUE(tree.get_root()->get_most_left_child()->has_childrens());
    EXPECT_FALSE(tree.get_root()->get_most_left_child()->is_leaf());
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_parent(), tree.get_root().get());
    EXPECT_EQ(tree.get_root()->get_most_left_child()->get_most_left_child(), tree.get_root()->get_most_left_child()->get_most_right_child());

    auto ab_child = tree.get_root()->get_most_left_child()->get_most_left_child();
    EXPECT_EQ(ab_child->get_data(), "AB");
    EXPECT_EQ(ab_child->get_type(), node_t::classic);
    EXPECT_TRUE(ab_child->has_childrens());
    EXPECT_FALSE(ab_child->is_leaf());
    EXPECT_EQ(ab_child->get_parent(), tree.get_root()->get_most_left_child().get());

    EXPECT_EQ(ab_child->get_most_left_child()->get_data(), "WGD1");
    EXPECT_EQ(ab_child->get_most_left_child()->get_type(), node_t::whole_duplication);
    EXPECT_EQ(ab_child->get_most_left_child()->get_parent(), ab_child.get());
    EXPECT_TRUE(ab_child->get_most_left_child()->has_childrens());
    EXPECT_EQ(ab_child->get_most_left_child()->get_most_left_child(), ab_child->get_most_left_child()->get_most_right_child());
    EXPECT_FALSE(ab_child->get_most_left_child()->is_leaf());

    EXPECT_EQ(ab_child->get_most_left_child()->get_most_left_child()->get_data(), "WGD2");
    EXPECT_EQ(ab_child->get_most_left_child()->get_most_left_child()->get_type(), node_t::whole_duplication);
    EXPECT_EQ(ab_child->get_most_left_child()->get_most_left_child()->get_parent(), ab_child->get_most_left_child().get());
    EXPECT_TRUE(ab_child->get_most_left_child()->get_most_left_child()->has_childrens());
    EXPECT_FALSE(ab_child->get_most_left_child()->get_most_left_child()->is_leaf());

    auto wgd2_child = ab_child->get_most_left_child()->get_most_left_child();
    EXPECT_EQ(wgd2_child->get_most_left_child()->get_data(), "WGD3");
    EXPECT_EQ(wgd2_child->get_most_left_child()->get_type(), node_t::whole_duplication);
    EXPECT_EQ(wgd2_child->get_most_left_child()->get_parent(), wgd2_child.get());
    EXPECT_TRUE(wgd2_child->get_most_left_child()->has_childrens());
    EXPECT_EQ(wgd2_child->get_most_left_child(), wgd2_child->get_most_right_child());
    EXPECT_TRUE(wgd2_child->get_most_left_child()->get_most_left_child()->is_leaf());
}

TEST_F(AlgorihmsTest, TwoGeneomesGetMulticolor) {
    //TEST CLASSIC TREE
    mcolor_t complete_color;
    mcolor_name_t mcolor_name;
    string_tree_t tree = string_tree_builder_t::build_tree("(A,B)AB");

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);

    EXPECT_EQ(mcolor_name.find(Mcolor(1))->second, "A");
    EXPECT_EQ(mcolor_name.find(Mcolor(2))->second, "B");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 2}))->second, "AB");
    EXPECT_EQ(Mcolor({1, 2}), complete_color);
    EXPECT_EQ(mcolor_name.size(), 3);

    //TEST TREE WITH BASIC wgd
    mcolor_name.clear();
    tree = string_tree_builder_t::build_tree("(A,B)AB");
    map_wgd_t wgds({all_wgds[0], all_wgds[3]});
    refine_tree_by_wgds(tree, wgds);

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);

    EXPECT_EQ(mcolor_name.find(Mcolor(1))->second, "A");
    EXPECT_EQ(mcolor_name.find(Mcolor(2))->second, "B");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1}))->second, "WGD1");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 2}))->second, "AB");
    EXPECT_EQ(Mcolor({1, 1, 2}), complete_color);
    EXPECT_EQ(mcolor_name.size(), 4);

    //TEST TREE WITH MULTIPLY WGDS on branch
    mcolor_name.clear();
    tree = string_tree_builder_t::build_tree("(A,B)AB");
    wgds = map_wgd_t({all_wgds[1]});
    refine_tree_by_wgds(tree, wgds);

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);

    EXPECT_EQ(mcolor_name.find(Mcolor(1))->second, "A");
    EXPECT_EQ(mcolor_name.find(Mcolor(2))->second, "B");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1}))->second, "WGD3");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1}))->second, "WGD2");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1, 1, 1, 1, 1}))->second, "WGD1");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1, 1, 1, 1, 1, 2}))->second, "AB");
    EXPECT_EQ(Mcolor({1, 1, 1, 1, 1, 1, 1, 1, 2}), complete_color);
    EXPECT_EQ(mcolor_name.size(), 6);
}

TEST_F(AlgorihmsTest, NGeneomesFullTreeGetMulticolor) {
    //TEST CLASSIC TREE
    mcolor_t complete_color;
    mcolor_name_t mcolor_name;
    string_tree_t tree = string_tree_builder_t::build_tree("(((A, B)AB, C) ABC, ((D, E)DE, F)DEF )");

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);

    for (auto const & info : genome_number) {
        EXPECT_EQ(mcolor_name.find(Mcolor(info.second))->second, info.first);
    }

    EXPECT_EQ(mcolor_name.find(Mcolor({1, 2}))->second, "AB");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 2, 3}))->second, "ABC");
    EXPECT_EQ(mcolor_name.find(Mcolor({4, 5}))->second, "DE");
    EXPECT_EQ(mcolor_name.find(Mcolor({4, 5, 6}))->second, "DEF");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 2, 3, 4, 5, 6}))->second, "");
    EXPECT_EQ(Mcolor({1, 2, 3, 4, 5, 6}), complete_color);
    EXPECT_EQ(mcolor_name.size(), 11);

    //TEST TREE WITH BASIC wgd
    mcolor_name.clear();
    tree = string_tree_builder_t::build_tree("(((A, B)AB, C) ABC, ((D, E)DE, F)DEF )");
    map_wgd_t wgds({all_wgds[0], all_wgds[3], all_wgds[4],all_wgds[5]});
    refine_tree_by_wgds(tree, wgds);

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);

    for (auto const & info : genome_number) {
        EXPECT_EQ(mcolor_name.find(Mcolor(info.second))->second, info.first);
    }

    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1}))->second, "WGD1");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 2}))->second, "AB");
    EXPECT_EQ(mcolor_name.find(Mcolor({3, 3}))->second, "WGD5");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1, 2, 2}))->second, "WGD6");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1, 2, 2, 3, 3}))->second, "ABC");
    EXPECT_EQ(mcolor_name.find(Mcolor({4, 5}))->second, "DE");
    EXPECT_EQ(mcolor_name.find(Mcolor({4, 5, 6}))->second, "DEF");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6}))->second, "");
    EXPECT_EQ(Mcolor({1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6}), complete_color);
    EXPECT_EQ(mcolor_name.size(), 14);
}

TEST_F(AlgorihmsTest, NGeneomesNonFullTreeGetMulticolor) {
    mcolor_t complete_color;
    genome_number.insert(std::make_pair("G", 7));

    //TEST CLASSIC TREE
    mcolor_name_t mcolor_name;
    string_tree_t tree = string_tree_builder_t::build_tree("(((A, B)AB, C, D)ABCD, (E,F)EF, G)");
    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);

    for (auto const & info : genome_number) {
        EXPECT_EQ(mcolor_name.find(Mcolor(info.second))->second, info.first);
    }

    EXPECT_EQ(mcolor_name.find(Mcolor({1, 2}))->second, "AB");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 2, 3, 4}))->second, "ABCD");
    EXPECT_EQ(mcolor_name.find(Mcolor({5, 6}))->second, "EF");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 2, 3, 4, 5, 6, 7}))->second, "");
    EXPECT_EQ(Mcolor({1, 2, 3, 4, 5, 6, 7}), complete_color);
    EXPECT_EQ(mcolor_name.size(), 11);

    //TEST TREE WITH BASIC wgd
    mcolor_name.clear();
    tree = string_tree_builder_t::build_tree("(((A, B)AB, C, D)ABCD, (E,F)EF, G)");
    map_wgd_t wgds({all_wgds[0], all_wgds[3], all_wgds[4], all_wgds[5], all_wgds[6], all_wgds[7]});
    refine_tree_by_wgds(tree, wgds);

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);

    for (auto const & info : genome_number) {
        EXPECT_EQ(mcolor_name.find(Mcolor(info.second))->second, info.first);
    }

    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1}))->second, "WGD1");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 2}))->second, "AB");
    EXPECT_EQ(mcolor_name.find(Mcolor({3, 3}))->second, "WGD7");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1, 2, 2}))->second, "WGD8");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1, 2, 2, 3, 3, 4}))->second, "ABCD");
    EXPECT_EQ(mcolor_name.find(Mcolor({5, 6}))->second, "EF");
    EXPECT_EQ(mcolor_name.find(Mcolor({1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7}))->second, "");
    EXPECT_EQ(Mcolor({1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7}), complete_color);
    EXPECT_EQ(mcolor_name.size(), 14);
}

TEST_F(AlgorihmsTest, TwoGeneomesGetMedians) {
    medians_color_t medians;
    mcolor_t complete_color;
    mcolor_name_t mcolor_name;
    name_mcolor_t name_mcolor;

    //TEST CLASSIC TREE
    string_tree_t tree = string_tree_builder_t::build_tree("(A,B)AB");

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);
    name_mcolor = swipe(mcolor_name);
    get_median_colors(tree, complete_color, name_mcolor, medians);
    EXPECT_EQ(medians.size(), 0);

    //TEST TREE WITH wgd
    medians.clear();
    mcolor_name.clear();
    map_wgd_t wgds({all_wgds[1]});
    tree = string_tree_builder_t::build_tree("(A,B)AB");
    refine_tree_by_wgds(tree, wgds);

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);
    name_mcolor = swipe(mcolor_name);
    get_median_colors(tree, complete_color, name_mcolor, medians);
    EXPECT_EQ(medians.size(), 3);
}

TEST_F(AlgorihmsTest, ThreeGeneomesGetMedians) {
    medians_color_t medians;
    mcolor_t complete_color;
    mcolor_name_t mcolor_name;
    name_mcolor_t name_mcolor;

    //TEST CLASSIC TREE
    string_tree_t tree = string_tree_builder_t::build_tree("((A,B)AB,C)");

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);
    name_mcolor = swipe(mcolor_name);
    get_median_colors(tree, complete_color, name_mcolor, medians);
    EXPECT_EQ(medians.size(), 1);
    EXPECT_EQ(std::get<0>(medians[0]), mcolor_t(1));
    EXPECT_EQ(std::get<1>(medians[0]), mcolor_t(2));
    EXPECT_EQ(std::get<2>(medians[0]), mcolor_t(3));

    //TEST TREE WITH wgd
    medians.clear();
    mcolor_name.clear();
    map_wgd_t wgds({all_wgds[0]});
    tree = string_tree_builder_t::build_tree("((A,B)AB,C)");
    refine_tree_by_wgds(tree, wgds);

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);
    name_mcolor = swipe(mcolor_name);
    get_median_colors(tree, complete_color, name_mcolor, medians);
    EXPECT_EQ(medians.size(), 2);
    EXPECT_EQ(std::get<0>(medians[0]), mcolor_t(1));
    EXPECT_EQ(std::get<1>(medians[0]), mcolor_t(1));
    EXPECT_EQ(std::get<2>(medians[0]), mcolor_t({2, 3}));
    EXPECT_EQ(std::get<0>(medians[1]), mcolor_t({1, 1}));
    EXPECT_EQ(std::get<1>(medians[1]), mcolor_t(2));
    EXPECT_EQ(std::get<2>(medians[1]), mcolor_t(3));
}

TEST_F(AlgorihmsTest, NGeneomesGetMedians) {
    medians_color_t medians;
    mcolor_t complete_color;
    mcolor_name_t mcolor_name;
    name_mcolor_t name_mcolor;
    genome_number.insert(std::make_pair("G", 7));

    //TEST CLASSIC TREE
    string_tree_t tree = string_tree_builder_t::build_tree("(((A, B)AB, C, D)ABCD, (E,F)EF, G)");
    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);
    name_mcolor = swipe(mcolor_name);
    get_median_colors(tree, complete_color, name_mcolor, medians);

    EXPECT_EQ(medians.size(), 2);
    EXPECT_EQ(std::get<0>(medians[0]), mcolor_t(1));
    EXPECT_EQ(std::get<1>(medians[0]), mcolor_t(2));
    EXPECT_EQ(std::get<2>(medians[0]), mcolor_t({3, 4, 5, 6, 7}));
    EXPECT_EQ(std::get<0>(medians[1]), mcolor_t(5));
    EXPECT_EQ(std::get<1>(medians[1]), mcolor_t(6));
    EXPECT_EQ(std::get<2>(medians[1]), mcolor_t({1, 2, 3, 4, 7}));

    //TEST TREE WITH wgd
    medians.clear();
    mcolor_name.clear();
    map_wgd_t wgds({all_wgds[0]});
    tree = string_tree_builder_t::build_tree("(((A, B)AB, C, D)ABCD, (E,F)EF, G)");
    refine_tree_by_wgds(tree, wgds);

    get_T_multicolors(tree, genome_number, mcolor_name, complete_color);
    name_mcolor = swipe(mcolor_name);
    get_median_colors(tree, complete_color, name_mcolor, medians);

    EXPECT_EQ(medians.size(), 3);
    EXPECT_EQ(std::get<0>(medians[0]), mcolor_t(1));
    EXPECT_EQ(std::get<1>(medians[0]), mcolor_t(1));
    EXPECT_EQ(std::get<2>(medians[0]), mcolor_t({2, 3, 4, 5, 6, 7}));
    EXPECT_EQ(std::get<0>(medians[1]), mcolor_t({1, 1}));
    EXPECT_EQ(std::get<1>(medians[1]), mcolor_t(2));
    EXPECT_EQ(std::get<2>(medians[1]), mcolor_t({3, 4, 5, 6, 7}));
    EXPECT_EQ(std::get<0>(medians[2]), mcolor_t(5));
    EXPECT_EQ(std::get<1>(medians[2]), mcolor_t(6));
    EXPECT_EQ(std::get<2>(medians[2]), mcolor_t({1, 1, 2, 3, 4, 7}));
}

//TODO write test on destroy on branches
