//
// Created by pavel on 1/12/16.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"
#include "graph/graph_pack/multicolors.hpp"

using namespace std;
using namespace structure;
using namespace structure::phyl_tree;
using namespace graph::graph_pack;

using mcolor_t = Mcolor;

using wgd_t = event::wgd<std::string>;
using string_tree_t = BinaryTree<std::string>;
using node_t = string_tree_t::node_t;
using node_ptr = string_tree_t::node_ptr;
using string_tree_builder_t = TreeBuilder<string_tree_t>;

using mcolor_name_t = std::map<mcolor_t, std::string>;
using name_mcolor_t = std::unordered_map<std::string, mcolor_t>;
using medians_color_t = std::vector<std::tuple<mcolor_t, mcolor_t, mcolor_t>>;

class MulticolorsTest : public testing::Test {
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


TEST_F(MulticolorsTest, TreeWithTwoGenomesNoWgdEmptyTarget) {
    string_tree_t tree = string_tree_builder_t::build_tree("(A,B)AB");
    Multicolors<mcolor_t> multicolors(tree, genome_number, "");

    EXPECT_EQ(multicolors.get_complete_color(), Mcolor({1, 2}));
    EXPECT_EQ(multicolors.get_root_color(), Mcolor(2));

    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor(1)));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor(2)));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor()));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 1})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 2})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 2, 3})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor(3)));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({4, 4, 5})));

    EXPECT_TRUE(multicolors.is_vec_T_consistent_color(Mcolor(1)));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor(2)));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor()));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 2})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 2, 3})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor(3)));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({4, 4, 5})));

    EXPECT_EQ(multicolors.get_median_colors().size(), 0);
    EXPECT_EQ(multicolors.get_complement_color(Mcolor({1, 2})), Mcolor());
    EXPECT_EQ(multicolors.get_complement_color(Mcolor(1)), Mcolor(2));
    EXPECT_EQ(multicolors.get_complement_color(Mcolor(2)), Mcolor(1));

    std::set<mcolor_t> tests({Mcolor(1), Mcolor(2)});
    std::set<mcolor_t> temp = multicolors.split_color_on_vtc_color(Mcolor(1));
    EXPECT_EQ(temp.size(), 1);
    EXPECT_TRUE(tests.count(*temp.begin()) != 0);

    temp = multicolors.split_color_on_vtc_color(Mcolor({1, 2}));
    EXPECT_EQ(temp.size(), 1);
    EXPECT_EQ(*temp.begin(), Mcolor({1, 2}));

    temp = multicolors.split_color_on_tc_color(Mcolor(1), 2);
    EXPECT_EQ(temp.size(), 1);
    EXPECT_TRUE(tests.count(*temp.begin()) != 0);

    EXPECT_EQ(multicolors.get_min_addit_color_for_tc(Mcolor(1)), Mcolor());
    EXPECT_EQ(multicolors.get_min_addit_color_for_tc(Mcolor({1, 2})), Mcolor({1, 2}));
}

TEST_F(MulticolorsTest, TreeWithTwoGenomesWgdEmptyTarget) {
    map_wgd_t wgds({all_wgds[1]});
    string_tree_t tree = string_tree_builder_t::build_tree("(A,B)AB");
    refine_tree_by_wgds(tree, wgds);

    Multicolors<mcolor_t> multicolors(tree, genome_number, "");

    EXPECT_EQ(multicolors.get_complete_color(), Mcolor({1, 1, 1, 1, 1, 1, 1, 1, 2}));
    EXPECT_EQ(multicolors.get_root_color(), Mcolor({1, 1, 1, 1, 1, 1, 1, 1}));

    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor(1)));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1, 1, 1, 1, 1, 2})));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({1, 1})));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1, 1, 1, 1, 2})));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1, 1})));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1, 1, 2})));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor(2)));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1, 1, 1, 1, 1, 1})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1, 1, 1, 1, 1, 1, 2})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1, 2, 2})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 1, 1, 1, 2, 2})));

    EXPECT_TRUE(multicolors.is_vec_T_consistent_color(Mcolor(1)));
    EXPECT_TRUE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1})));
    EXPECT_TRUE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1, 1})));
    EXPECT_TRUE(multicolors.is_vec_T_consistent_color(Mcolor(2)));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1, 1, 1, 1, 1, 1})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1, 1, 1, 1, 1, 2})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1, 1, 1, 1, 2})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1, 1, 1, 1, 1, 1, 2})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1, 1, 2})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1, 2, 2})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1, 1, 1, 2, 2})));

    std::set<mcolor_t> tests({Mcolor({1, 1, 1, 1}), Mcolor({1, 1}), Mcolor({1}), Mcolor({2}), Mcolor({1, 1, 1, 1, 2})});
    std::set<mcolor_t> temp = multicolors.split_color_on_vtc_color(Mcolor(1));
    EXPECT_EQ(temp.size(), 1);
    EXPECT_TRUE(tests.count(*temp.begin()) != 0);

    temp = multicolors.split_color_on_vtc_color(Mcolor({1, 1, 1}));
    EXPECT_EQ(temp.size(), 2);
    EXPECT_TRUE(tests.count(*temp.begin()) != 0);
    EXPECT_TRUE(tests.count(*temp.rbegin()) != 0);

    temp = multicolors.split_color_on_vtc_color(Mcolor({1, 1, 1, 1, 1, 1, 1}));
    EXPECT_EQ(temp.size(), 3);
    for (auto const & color : temp) {
        EXPECT_TRUE(tests.count(color) != 0);
    }

    temp = multicolors.split_color_on_vtc_color(Mcolor({1, 1, 1, 1, 1, 1, 1, 2}));
    EXPECT_EQ(temp.size(), 4);
    for (auto const & color : temp) {
        EXPECT_TRUE(tests.count(color) != 0);
    }

    temp = multicolors.split_color_on_tc_color(Mcolor({1, 1, 1, 1, 1, 1, 1, 2}));
    EXPECT_EQ(temp.size(), 1);
    EXPECT_EQ(*temp.begin(), Mcolor({1, 1, 1, 1, 1, 1, 1, 2}));

    temp = multicolors.split_color_on_tc_color(Mcolor({1, 1, 1}));
    EXPECT_EQ(temp.size(), 1);

    temp = multicolors.split_color_on_tc_color(Mcolor({1, 1, 1}), 2);
    EXPECT_EQ(temp.size(), 2);
    EXPECT_TRUE(tests.count(*temp.begin()) != 0);
    EXPECT_TRUE(tests.count(*temp.rbegin()) != 0);

    temp = multicolors.split_color_on_tc_color(Mcolor({1, 1, 1, 1, 1, 1, 1}), 3);
    EXPECT_EQ(temp.size(), 3);
    for (auto const & color : temp) {
        EXPECT_TRUE(tests.count(color) != 0);
    }

    temp = multicolors.split_color_on_tc_color(Mcolor({1, 1, 1, 1, 1, 2}), 3);
    EXPECT_EQ(temp.size(), 2);
    for (auto const & color : temp) {
        EXPECT_TRUE(tests.count(color) != 0);
    }

    temp = multicolors.split_color_on_tc_color(Mcolor({1, 1, 1, 1, 1, 2}), 2);
    EXPECT_EQ(temp.size(), 2);
    EXPECT_TRUE(tests.count(*temp.rbegin()) != 0);
}

TEST_F(MulticolorsTest, TreeWithThreeGenomesNoWgdEmptyTarget) {
    string_tree_t tree = string_tree_builder_t::build_tree("((A,B)AB, C)");
    Multicolors<mcolor_t> multicolors(tree, genome_number, "");

    EXPECT_EQ(multicolors.get_complete_color(), Mcolor({1, 2, 3}));
    EXPECT_EQ(multicolors.get_root_color(), Mcolor({1, 2}));

    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor(1)));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor(2)));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor(3)));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({1, 2})));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({2, 3})));
    EXPECT_TRUE(multicolors.is_T_consistent_color(Mcolor({1, 3})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor()));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 1})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 1, 2})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({1, 2, 3})));
    EXPECT_FALSE(multicolors.is_T_consistent_color(Mcolor({4, 4, 5})));

    EXPECT_TRUE(multicolors.is_vec_T_consistent_color(Mcolor(1)));
    EXPECT_TRUE(multicolors.is_vec_T_consistent_color(Mcolor(2)));
    EXPECT_TRUE(multicolors.is_vec_T_consistent_color(Mcolor(3)));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 2})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor()));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({2, 3})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 3})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 1})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({1, 2, 3})));
    EXPECT_FALSE(multicolors.is_vec_T_consistent_color(Mcolor({4, 4, 5})));
    //todo continue here
}


/*TEST_F(MulticolorsTest, TreeWithThreeGenomesNoWgdEmptyTarget) {
    //todo continue here
}*/

/*TEST_F(MulticolorsTest, TreeWithThreeGenomesNoWgdTarget) {
    //todo continue here
}*/

TEST_F(MulticolorsTest, TreeWithThreeGenomesWgdEmptyTarget) {
    //todo continue here
}

TEST_F(MulticolorsTest, TreeWithThreeGenomesWgdTarget) {
    //todo continue here
}

TEST_F(MulticolorsTest, FullTreeWithNGenomesNoWgdEmptyTarget) {
    //todo continue here
}

TEST_F(MulticolorsTest, FullTreeWithNGenomesNoWgdTarget) {
    //todo continue here
}

TEST_F(MulticolorsTest, FullTreeWithNGenomesWgdEmptyTarget) {
    //todo continue here
}

TEST_F(MulticolorsTest, FullTreeWithNGenomesWgdTarget) {
    //todo continue here
}
