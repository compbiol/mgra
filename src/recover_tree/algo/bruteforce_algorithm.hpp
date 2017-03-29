#ifndef BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP_
#define BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP_

#include "recover_tree_algorithm.hpp"
#include "statistics_based_recover_tree_algorithm.hpp"

#include "structures/phyl_tree/branch.hpp"
#include "graph/graph_pack.hpp"

namespace algo {

template<class graph_pack_t>
struct BruteforceRecoverTreeAlgorithm :
        StatisticsBasedRecoverTreeAlgorithm<graph_pack_t> {
    using statistic_producer_ptr = std::shared_ptr<StatisticsProducer<graph_pack_t> > const &;

    using mcolor_t = typename RecoverTreeAlgorithm<graph_pack_t>::mcolor_t;
    using tree_t = typename RecoverTreeAlgorithm<graph_pack_t>::tree_t;
    using tree_ptr = typename RecoverTreeAlgorithm<graph_pack_t>::tree_ptr;
    using node_t = typename tree_t::node_t;
    using node_ptr = typename tree_t::node_ptr;

    using branch_t = structure::phyl_tree::Branch<mcolor_t>;

    using class_t = std::pair<std::vector<branch_t>, size_t>;
    using class_vector = std::vector<class_t>;
    using statistic_t = std::pair<branch_t, size_t>;


    BruteforceRecoverTreeAlgorithm(graph_pack_t &graph_pack,
                                   statistic_producer_ptr const & statistic_producer) :
            StatisticsBasedRecoverTreeAlgorithm<graph_pack_t>(graph_pack, statistic_producer) {
    }

    std::vector<tree_ptr> recover_trees() {
        auto statistics = this->m_statistic_producer->make_statistics();
        auto classes = get_branch_sets(statistics);
        return build_trees(classes, statistics);
    }

private:
    /**
    * Checks if the branch conflicts the class of branches, adds it if not
    * @return true if the branch was added
    */
    bool screen_branch(class_t &cls, statistic_t const & statistic) {
        auto const &new_branch = statistic.first;
        bool new_branch_conflicts = false;
        for (auto &branch: cls.first) {
            if (branch.do_intersect(new_branch)) {
                new_branch_conflicts = true;
                break;
            }
        }
        if (!new_branch_conflicts) {
            cls.first.push_back(new_branch);
            cls.second += statistic.second;
            return true;
        }
        return false;
    }

    class_vector get_branch_sets(std::vector<statistic_t> & branch_statistics) {
        class_vector tree_classes;

        for (size_t i = 0; i != branch_statistics.size(); ++i) {
            auto const &statistic = branch_statistics[i];

            bool has_been_added = false;
            for (auto &cls: tree_classes) {
                bool screening_result = screen_branch(cls, statistic);
                has_been_added = has_been_added || screening_result;
            }

            if (!has_been_added) {
                tree_classes.push_back(class_t({statistic.first}, statistic.second));
                class_t &new_class = tree_classes.back();

                // So we won't forget to add the branches already screened, but consistent with new class
                for (size_t j = 0; j != i; ++j) {
                    screen_branch(new_class, branch_statistics[j]);
                }
            }
        }

        return tree_classes;
    }

    std::vector<tree_ptr> build_trees(std::vector<class_t> & tree_classes, std::vector<statistic_t> const & branch_statistics) {
        using namespace structure::phyl_tree;

        // Descending by class score
        std::sort(std::begin(tree_classes), std::end(tree_classes), [](class_t const &left, class_t const &right) {
            return left.second > right.second;
        });

        // Collect the colors appearing in the tree to speed up the strictness checking
        std::set<mcolor_t> appearing_colors;
        for (auto const &statistic: branch_statistics) {
            auto left_color = statistic.first.left;
            appearing_colors.insert(left_color);
        }

        std::vector<tree_ptr> results;
        results.reserve(tree_classes.size());
        /*std::transform(std::begin(tree_classes), std::end(tree_classes), std::back_inserter(results),
                       [&appearing_colors](class_t &cls_to_fold) {
                           //auto root_node = TreeBuilder<tree_t>::build_tree(cls_to_fold.first, cfg::get().priority_name.size());
                           ;//does_need_to_prune_node(root_node, appearing_colors); FIXME ADD PRUNING
                           return nullptr; //std::make_shared<tree_t>(root_node); //TODO come back
                       });*/
        return results;
    }

};

}

#endif