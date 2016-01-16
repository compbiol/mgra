//
// Created by pavel on 1/14/16.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"

using namespace structure;
using namespace structure::phyl_tree;

using string_tree_t = BinaryTree<std::string>;
using string_tree_builder_t = TreeBuilder<string_tree_t>;
using string_node_t = string_tree_t::node_t;
using string_node_ptr = string_tree_t::node_ptr;
using string_printer = NewickTreePrinter<string_tree_t>;

using wgd_t = event::wgd<std::string>;
using map_wgd_t = std::map<std::pair<std::string, std::string>, wgd_t>;

class NewickTreePrinterTest : public testing::Test {
    virtual void SetUp() override {
        all_wgds.insert(std::make_pair(std::make_pair("AB", "A"), wgd_t("AB", "A", std::vector<std::string>({"WGD1"}))));
        all_wgds.insert(std::make_pair(std::make_pair("root", "AB"), wgd_t("root", "AB", std::vector<std::string>({"WGD2"}))));
        all_wgds.insert(std::make_pair(std::make_pair("AB", "D"), wgd_t("AB", "D", std::vector<std::string>({"WGD3"}))));
        all_wgds.insert(std::make_pair(std::make_pair("DR", "D"), wgd_t("DR", "D", std::vector<std::string>({"WGD4", "WGD5"}))));
        all_wgds.insert(std::make_pair(std::make_pair("DRC", "DR"), wgd_t("DRC", "DR", std::vector<std::string>({"WGD6"}))));
        all_wgds.insert(std::make_pair(std::make_pair("ABC", "B"), wgd_t("ABC", "B", std::vector<std::string>({"WGD7", "WGD8", "WGD9"}))));
        all_wgds.insert(std::make_pair(std::make_pair("ABC", "C"), wgd_t("ABC", "C", std::vector<std::string>({"WGD10"}))));
        all_wgds.insert(std::make_pair(std::make_pair("ABC", "AB"), wgd_t("ABC", "AB", std::vector<std::string>({"WGD11"}))));
        all_wgds.insert(std::make_pair(std::make_pair("FDEGHN", "GHN"), wgd_t("FDEGHN", "GHN", std::vector<std::string>({"WGD12"}))));
        all_wgds.insert(std::make_pair(std::make_pair("FDEGHN", "FDE"), wgd_t("FDEGHN", "FDE", std::vector<std::string>({"WGD13"}))));
        all_wgds.insert(std::make_pair(std::make_pair("ABCDEF", "AB"), wgd_t("ABCDEF", "AB", std::vector<std::string>({"WGD14"}))));
        all_wgds.insert(std::make_pair(std::make_pair("ABCDEF", "CD"), wgd_t("ABCDEF", "CD", std::vector<std::string>({"WGD15"}))));
        all_wgds.insert(std::make_pair(std::make_pair("ABCDEF", "EF"), wgd_t("ABCDEF", "EF", std::vector<std::string>({"WGD16"}))));
    }

protected:
    map_wgd_t all_wgds;
};

TEST_F(NewickTreePrinterTest, BasicTreeTwoGenomes) {
    std::vector<std::string> in_strs({"(A, B) ", "(A, B)AB:3", "(Dog, Rat)", "(Dog, Rat):3 "});
    std::vector<std::string> out_strs({"(A, B);", "(A, B)AB;", "(Dog, Rat);", "(Dog, Rat);"});

    assert(in_strs.size() == out_strs.size());
    for (size_t i = 0; i < in_strs.size(); ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(in_strs[i]);
        std::stringstream out;
        string_printer newick_printer(out);
        newick_printer.print_tree(tree);
        EXPECT_EQ(out.str(), out_strs[i]);
    }
}

TEST_F(NewickTreePrinterTest, BasicTreeTwoGenomesWithWGD) {
    std::vector<std::string> in_strs({"(A, B)AB:3", " ( A , B ) AB : 3 "});
    std::vector<std::string> out_strs({"(A, B)AB;", "(A, B)AB;"});

    assert(in_strs.size() == out_strs.size());
    for (size_t i = 0; i < in_strs.size(); ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(in_strs[i]);
        refine_tree_by_wgds(tree, all_wgds);
        std::stringstream out;
        string_printer newick_printer(out);
        newick_printer.print_tree(tree);
        EXPECT_EQ(out.str(), out_strs[i]);
    }

}

TEST_F(NewickTreePrinterTest, BasicTreeThreeGenomes) {
    std::vector<std::string> in_strs(
            {"((A, B), C)", "((A, B) : 3, C)", "((A, B), C) : 4", "((A, B) : 3, C) : 5", "((Dog, Rat)DR, C)",
             "((Dog, Rat), C)DRC", "((Dog, Rat)DR, C)DRC", "((Dog, Rat)DR:3, C)DRC:4"});
    std::vector<std::string> out_strs(
            {"((A, B), C);", "((A, B), C);", "((A, B), C);", "((A, B), C);", "((Dog, Rat)DR, C);",
             "((Dog, Rat), C)DRC;", "((Dog, Rat)DR, C)DRC;", "((Dog, Rat)DR, C)DRC;"});

    assert(in_strs.size() == out_strs.size());
    for (size_t i = 0; i < in_strs.size(); ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(in_strs[i]);
        std::stringstream out;
        string_printer newick_printer(out);
        newick_printer.print_tree(tree);
        EXPECT_EQ(out.str(), out_strs[i]);
    }
}

TEST_F(NewickTreePrinterTest, BasicTreeThreeGenomesWithWGD) {
    std::vector<std::string> in_strs(
            {"((A, B), C)", "((A, B) : 3, C)", "((A, B), C) : 4", "((A, B) : 3, C) : 5", "((Dog, Rat)DR, C)",
             "((Dog, Rat), C)DRC", "((Dog, Rat)DR, C)DRC", "((Dog, Rat)DR:3, C)DRC:4"});
    std::vector<std::string> out_strs(
            {"((A, B), C);", "((A, B), C);", "((A, B), C);", "((A, B), C);", "((Dog, Rat)DR, C);",
             "((Dog, Rat), C)DRC;", "((Dog, Rat)DR, C)DRC;", "((Dog, Rat)DR, C)DRC;"});

    assert(in_strs.size() == out_strs.size());
    for (size_t i = 0; i < in_strs.size(); ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(in_strs[i]);
        refine_tree_by_wgds(tree, all_wgds);
        std::stringstream out;
        string_printer newick_printer(out);
        newick_printer.print_tree(tree);
        EXPECT_EQ(out.str(), out_strs[i]);
    }
}

TEST_F(NewickTreePrinterTest, GeneralTreeThreeGenomes) {
    std::vector<std::string> in_strs({"(A, B, C)", "(A, B, C) : 4", "(A, B, C)ABC", "(A, B, C)ABC:4"});
    std::vector<std::string> out_strs({"(A, B, C);", "(A, B, C);", "(A, B, C)ABC;", "(A, B, C)ABC;"});

    assert(in_strs.size() == out_strs.size());
    for (size_t i = 0; i < in_strs.size(); ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(in_strs[i]);
        std::stringstream out;
        string_printer newick_printer(out);
        newick_printer.print_tree(tree);
        EXPECT_EQ(out.str(), out_strs[i]);
    }
}

TEST_F(NewickTreePrinterTest, GeneralTreeThreeGenomesWithWGD) {
    std::vector<std::string> in_strs({"(A, B, C)", "(A, B, C) : 4", "(A, B, C)ABC", "(A, B, C)ABC:4"});
    std::vector<std::string> out_strs({"(A, B, C);", "(A, B, C);", "(A, B, C)ABC;", "(A, B, C)ABC;"});

    assert(in_strs.size() == out_strs.size());
    for (size_t i = 0; i < in_strs.size(); ++i) {
        string_tree_t tree = string_tree_builder_t::build_tree(in_strs[i]);
        refine_tree_by_wgds(tree, all_wgds);
        std::stringstream out;
        string_printer newick_printer(out);
        newick_printer.print_tree(tree);
        EXPECT_EQ(out.str(), out_strs[i]);
    }
}

TEST_F(NewickTreePrinterTest, BasicTreeNGenomes) {
    string_tree_t tree = string_tree_builder_t::build_tree(
            "(((A,B)AB,C)ABC,((F,(D,E) : 4)FDE,(G,(H,N)HN)GHN:3)FDEGHN)");
    std::stringstream out;
    string_printer newick_printer(out);
    newick_printer.print_tree(tree);
    EXPECT_EQ(out.str(), "(((A, B)AB, C)ABC, ((F, (D, E))FDE, (G, (H, N)HN)GHN)FDEGHN);");
}

TEST_F(NewickTreePrinterTest, BasicTreeNGenomesWithWGD) {
    string_tree_t tree = string_tree_builder_t::build_tree(
            "(((A,B)AB,C)ABC,((F,(D,E) : 4)FDE,(G,(H,N)HN)GHN:3)FDEGHN)");
    refine_tree_by_wgds(tree, all_wgds);
    std::stringstream out;
    string_printer newick_printer(out);
    newick_printer.print_tree(tree);
    EXPECT_EQ(out.str(), "(((A, B)AB, C)ABC, ((F, (D, E))FDE, (G, (H, N)HN)GHN)FDEGHN);");
}

TEST_F(NewickTreePrinterTest, GeneralTreeNGenomes) {
    string_tree_t tree = string_tree_builder_t::build_tree(
            "((((A,B)AB, (C, D)CD:3,(E,F):4)ABCDEF,G), (H,(M,N,R)MNR)MNRH:10,(O,P):4)");
    std::stringstream out;
    string_printer newick_printer(out);
    newick_printer.print_tree(tree);
    EXPECT_EQ(out.str(), "((((A, B)AB, (C, D)CD, (E, F))ABCDEF, G), (H, (M, N, R)MNR)MNRH, (O, P));");
}

TEST_F(NewickTreePrinterTest, GeneralTreeNGenomesWithWGD) {
    string_tree_t tree = string_tree_builder_t::build_tree(
            "((((A,B)AB, (C, D)CD:3,(E,F):4)ABCDEF,G), (H,(M,N,R)MNR)MNRH:10,(O,P):4)");
    refine_tree_by_wgds(tree, all_wgds);
    std::stringstream out;
    string_printer newick_printer(out);
    newick_printer.print_tree(tree);
    EXPECT_EQ(out.str(), "((((A, B)AB, (C, D)CD, (E, F))ABCDEF, G), (H, (M, N, R)MNR)MNRH, (O, P));");
}

