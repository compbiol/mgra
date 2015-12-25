//
// Created by pavel on 12/21/15.
//

#include "gtest/gtest.h"
#include "defined_1.hpp"
#include "structures/mcolor.hpp"

using namespace structure;

class McolorTest : public testing::Test {
    virtual void SetUp() override {
        colors.push_back(Mcolor());                                                 //0

        colors.push_back(Mcolor(1)); colors[1].insert(1); colors[1].insert(1);      //1
        colors[1].insert(2); colors[1].insert(2); colors[1].insert(3);
        colors[1].insert(3);
        colors.push_back(Mcolor(1, 3)); colors[2].insert(2, 2);                     //2
        colors[2].insert(3, 2);
        colors.push_back(Mcolor()); colors[3].insert(1, 3); colors[3].insert(2, 2); //3
        colors[3].insert(3, 2);
        colors.push_back(Mcolor(1)); colors[4].insert(2); colors[4].insert(3);      //4

        colors.push_back(Mcolor({1, 2, 3}));                                        //5
        colors.push_back(Mcolor({4, 5, 6}));                                        //6
        colors.push_back(Mcolor({1, 2, 3, 4, 5, 6}));                               //7
        colors.push_back(Mcolor({1, 1, 2, 2, 3, 3}));                               //8
        colors.push_back(Mcolor({4, 4, 5, 5, 6, 6}));                               //9
        colors.push_back(Mcolor({1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6}));             //10

        colors.push_back(Mcolor({1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6})); //11
        colors.push_back(Mcolor({1, 2, 3, 4, 4, 5, 5, 6, 6}));                            //12
        colors.push_back(Mcolor({1, 1, 1, 7, 7, 8}));                                     //13
        colors.push_back(Mcolor({1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6})); //14
    }

protected:
    std::vector<Mcolor> colors;
};

TEST_F(McolorTest, BasicConstructors) {
    EXPECT_EQ(Mcolor(), Mcolor({}));
    EXPECT_EQ(Mcolor(1), Mcolor({1}));
    EXPECT_EQ(Mcolor(1, 2), Mcolor({1, 1}));
}

TEST_F(McolorTest, ConstructorByMcolorsIntersect) {
    //intersection of two disjoint simply sets
    Mcolor intersect(colors[5], colors[6], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[0]);

    intersect = Mcolor(colors[6], colors[5], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[0]);

    //intersection of two disjoint general sets
    intersect = Mcolor(colors[8], colors[9], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[0]);

    //intersection of one simple set with the same set
    intersect = Mcolor(colors[5], colors[5], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[5]);

    //intersection of one general set with the same set
    intersect = Mcolor(colors[1], colors[1], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[1]);

    //intersection of one simple set with the empty set
    intersect = Mcolor(colors[5], colors[0], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[0]);

    //intersection of one general set with the empty set
    intersect = Mcolor(colors[1], colors[0], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[0]);

    //intersection of two simply sets (full includes)
    intersect = Mcolor(colors[5], colors[7], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[5]);

    //intersection of two general sets (full includes)
    intersect = Mcolor(colors[8], colors[10], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[8]);

    //intersection of two simply sets
    intersect = Mcolor(Mcolor({4, 8, 9}), colors[6], Mcolor::Intersection);
    EXPECT_EQ(intersect, Mcolor(4));

    //intersection of two general sets
    intersect = Mcolor(colors[13], colors[8], Mcolor::Intersection);
    EXPECT_EQ(intersect, Mcolor(1, 2));

    //intersection of a general set and a simply set
    intersect = Mcolor(colors[7], colors[10], Mcolor::Intersection);
    EXPECT_EQ(intersect, colors[7]);
}

TEST_F(McolorTest, ConstructorByMcolorsUnion) {
    //union of two disjoint simply sets
    Mcolor result(colors[5], colors[6], Mcolor::Union);
    EXPECT_EQ(result, colors[7]);

    //union of two disjoint general sets
    result = Mcolor(colors[8], colors[9], Mcolor::Union);
    EXPECT_EQ(result, colors[10]);

    result = Mcolor(colors[9], colors[8], Mcolor::Union);
    EXPECT_EQ(result, colors[10]);

    //union of one simple set with the same set
    result = Mcolor(colors[5], colors[5], Mcolor::Union);
    EXPECT_EQ(result, colors[8]);

    result = Mcolor(colors[6], colors[6], Mcolor::Union);
    EXPECT_EQ(result, colors[9]);

    //union of one general set with the same set
    result = Mcolor(colors[1], colors[1], Mcolor::Union);
    EXPECT_EQ(result, Mcolor({1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3}));

    //union of one simple set with the empty set
    result = Mcolor(colors[5], colors[0], Mcolor::Union);
    EXPECT_EQ(result, colors[5]);

    //union of one general set with the empty set
    result = Mcolor(colors[1], colors[0], Mcolor::Union);
    EXPECT_EQ(result, colors[1]);

    //union of two simply sets (full includes)
    result = Mcolor(colors[6], colors[7], Mcolor::Union);
    EXPECT_EQ(result, colors[12]);

    //union of two general sets (full includes)
    result = Mcolor(colors[10], colors[8], Mcolor::Union);
    EXPECT_EQ(result, colors[11]);

    //union of two simply sets
    result = Mcolor(Mcolor({4, 8, 9}), colors[6], Mcolor::Union);
    EXPECT_EQ(result, Mcolor({4, 4, 5, 6, 8, 9}));

    //union of two general sets
    result = Mcolor(Mcolor({1, 1, 1, 7, 7, 8}), colors[8], Mcolor::Union);
    EXPECT_EQ(result, Mcolor({1, 1, 1, 1, 1, 2, 2, 3, 3, 7, 7, 8}));

    //union of a general set and a simply set
    result = Mcolor(colors[7], colors[10], Mcolor::Union);
    EXPECT_EQ(result, colors[14]);

    result = Mcolor(colors[10], colors[7], Mcolor::Union);
    EXPECT_EQ(result, colors[14]);
}

TEST_F(McolorTest, ConstructorByMcolorsDifferences) {
    Mcolor result;

    //difference of two disjoint simply sets
    result = Mcolor(colors[5], colors[6], Mcolor::Difference);
    EXPECT_EQ(result, colors[5]);

    result = Mcolor(colors[6], colors[5], Mcolor::Difference);
    EXPECT_EQ(result, colors[6]);

    //difference of two disjoint general sets
    result = Mcolor(colors[8], colors[9], Mcolor::Difference);
    EXPECT_EQ(result, colors[8]);

    result = Mcolor(colors[9], colors[8], Mcolor::Difference);
    EXPECT_EQ(result, colors[9]);

    //difference of one simple set with the same set
    result = Mcolor(colors[5], colors[5], Mcolor::Difference);
    EXPECT_EQ(result, colors[0]);

    //difference of one general set with the same set
    result = Mcolor(colors[1], colors[1], Mcolor::Difference);
    EXPECT_EQ(result, colors[0]);

    result = Mcolor(colors[8], colors[8], Mcolor::Difference);
    EXPECT_EQ(result, colors[0]);

    //difference of one simple set with the empty set
    result = Mcolor(colors[5], colors[0], Mcolor::Difference);
    EXPECT_EQ(result, colors[5]);

    result = Mcolor(colors[0], colors[5], Mcolor::Difference);
    EXPECT_EQ(result, colors[0]);

    //difference of one general set with the empty set
    result = Mcolor(colors[1], colors[0], Mcolor::Difference);
    EXPECT_EQ(result, colors[1]);

    result = Mcolor(colors[0], colors[1], Mcolor::Difference);
    EXPECT_EQ(result, colors[0]);

    //difference of two simply sets (full includes)
    result = Mcolor(colors[5], colors[7], Mcolor::Difference);
    EXPECT_EQ(result, colors[0]);

    result = Mcolor(colors[7], colors[5], Mcolor::Difference);
    EXPECT_EQ(result, colors[6]);

    //difference of two general sets (full includes)
    result = Mcolor(colors[9], colors[10], Mcolor::Difference);
    EXPECT_EQ(result, colors[0]);

    result = Mcolor(colors[10], colors[9], Mcolor::Difference);
    EXPECT_EQ(result, colors[8]);

    //difference of two simply sets
    result = Mcolor(Mcolor({4, 8, 9}), colors[6], Mcolor::Difference);
    EXPECT_EQ(result, Mcolor({8, 9}));

    result = Mcolor(colors[6], Mcolor({4, 8, 9}), Mcolor::Difference);
    EXPECT_EQ(result, Mcolor({5, 6}));

    //difference of two general sets
    result = Mcolor(Mcolor({1, 1, 1, 7, 7, 8}), colors[8], Mcolor::Difference);
    EXPECT_EQ(result, Mcolor({1, 7, 7, 8}));

    result = Mcolor(colors[8], Mcolor({1, 1, 1, 7, 7, 8}), Mcolor::Difference);
    EXPECT_EQ(result, Mcolor({2, 2, 3, 3}));

    //difference of a general set and a simply set
    result = Mcolor(colors[10], colors[7], Mcolor::Difference);
    EXPECT_EQ(result, colors[7]);

    result = Mcolor(colors[7], colors[10], Mcolor::Difference);
    EXPECT_EQ(result, colors[0]);
}

TEST_F(McolorTest, InsertFunc) {
    EXPECT_EQ(colors[1], colors[2]);
    EXPECT_EQ(colors[2], colors[3]);
    EXPECT_EQ(colors[1], colors[3]);
    EXPECT_EQ(colors[4], colors[5]);
}

TEST_F(McolorTest, ClearFunc) {
    colors[0].clear();
    EXPECT_EQ(colors[0], colors[0]);

    colors[1].clear();
    EXPECT_EQ(colors[1], colors[0]);

    colors[7].clear();
    EXPECT_EQ(colors[7], colors[0]);

    colors[10].clear();
    EXPECT_EQ(colors[10], colors[0]);
}

TEST_F(McolorTest, OneToOneFunc) {
    EXPECT_TRUE(colors[0].is_one_to_one_match());
    EXPECT_FALSE(colors[1].is_one_to_one_match());
    EXPECT_TRUE(colors[7].is_one_to_one_match());
    EXPECT_FALSE(colors[10].is_one_to_one_match());

    EXPECT_FALSE(Mcolor({1, 2, 3, 3}).is_one_to_one_match());
    EXPECT_FALSE(Mcolor({1, 2, 2, 3}).is_one_to_one_match());
}

TEST_F(McolorTest, IncludesFunc) {
    EXPECT_TRUE(colors[0].includes(colors[0]));
    EXPECT_TRUE(colors[7].includes(colors[7]));
    EXPECT_TRUE(colors[10].includes(colors[10]));

    EXPECT_TRUE(colors[7].includes(colors[0]));
    EXPECT_FALSE(colors[0].includes(colors[7]));

    EXPECT_TRUE(colors[10].includes(colors[0]));
    EXPECT_FALSE(colors[0].includes(colors[10]));

    EXPECT_TRUE(colors[7].includes(colors[5]));
    EXPECT_FALSE(colors[5].includes(colors[7]));

    EXPECT_TRUE(colors[7].includes(colors[6]));
    EXPECT_FALSE(colors[6].includes(colors[7]));

    EXPECT_TRUE(colors[10].includes(colors[8]));
    EXPECT_FALSE(colors[8].includes(colors[10]));

    EXPECT_TRUE(colors[10].includes(colors[9]));
    EXPECT_FALSE(colors[9].includes(colors[10]));

    EXPECT_TRUE(colors[10].includes(colors[7]));
    EXPECT_FALSE(colors[7].includes(colors[10]));

    EXPECT_TRUE(colors[10].includes(Mcolor{1, 2, 2, 3}));
}

TEST_F(McolorTest, DefinedFunc) {
    EXPECT_FALSE(colors[0].defined(1));
    EXPECT_FALSE(colors[0].defined(5));
    EXPECT_FALSE(colors[0].defined(7));

    EXPECT_TRUE(colors[1].defined(1));
    EXPECT_FALSE(colors[1].defined(5));
    EXPECT_FALSE(colors[1].defined(7));

    EXPECT_TRUE(colors[7].defined(1));
    EXPECT_TRUE(colors[7].defined(4));
    EXPECT_FALSE(colors[7].defined(7));
}

TEST_F(McolorTest, EmptyFunc) {
    EXPECT_TRUE(colors[0].empty());
    EXPECT_FALSE(colors[1].empty());
    EXPECT_FALSE(colors[5].empty());
    EXPECT_FALSE(colors[7].empty());
    EXPECT_FALSE(colors[10].empty());
}

TEST_F(McolorTest, SizeFunc) {
    EXPECT_EQ(colors[0].size(), 0);
    EXPECT_EQ(colors[1].size(), 7);
    EXPECT_EQ(colors[7].size(), 6);
    EXPECT_EQ(colors[10].size(), 12);
}