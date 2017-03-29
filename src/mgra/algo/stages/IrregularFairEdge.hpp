//
// Created by pavel on 1/15/16.
//

#ifndef MGRA_IRREGULAR_FAIR_EDGE_HPP
#define MGRA_IRREGULAR_FAIR_EDGE_HPP

//todo go to fair edge (Now only for WGD for two genomes)

namespace algo {

template<class graph_pack_t>
struct ProcessIrregularFairEdge : public algo::AbsStage<graph_pack_t> {

    using mcolor_t = typename graph_pack_t::mcolor_t;

    using edge_t = typename graph_pack_t::edge_t;
    using arc_t = typename graph_pack_t::arc_t;
    using mularcs_t = typename graph_pack_t::mularcs_t;
    using twobreak_t = typename graph_pack_t::twobreak_t;

    explicit ProcessIrregularFairEdge(size_t max_round)
            : algo::AbsStage<graph_pack_t>("Process irregular fair edges", "irregular_fair_edge", round_stage_t,
                                           max_round) {
    }

    bool run(graph_pack_t &graph_pack) override;


private:
    DECL_LOGGER("ProcessIrregularFairEdge");
};

template<class graph_pack_t>
bool ProcessIrregularFairEdge<graph_pack_t>::run(graph_pack_t &graph_pack) {
    bool isChanged = false;
    size_t number_rear = 0; // number of rearrangements

    do {
        number_rear = 0;
        //TODO modify mularcs for get on O(1) interesting vertex

        for (vertex_t const &x : graph_pack.graph) {
            mularcs_t const &mularcs = graph_pack.get_all_adjacent_multiedges(x);

            bool found = false;
            for (auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
                vertex_t const &y = im->first; // Q == im->second - color of central edge

                if (y == Infty) {
                    if (graph_pack.multicolors.is_T_consistent_color(im->second) &&
                        !graph_pack.multicolors.is_vec_T_consistent_color(im->second)) {
                        mularcs_t mularcs_x = graph_pack.get_all_adjacent_multiedges(x);
                        mularcs_x.erase(y);
                        for (auto it = mularcs_x.cbegin(); (it != mularcs_x.cend()); ++it) {
                            graph_pack.apply(twobreak_t(x, it->first, Infty, Infty, it->second));
                            found = true;
                            ++number_rear;
                        }
                    }
                }
            }
        }

        if (number_rear != 0) {
            isChanged = true;
        }
    } while (number_rear > 0);

    return isChanged;
}

}
#endif //MGRA_IRREGULAR_FAIR_EDGE_HPP
