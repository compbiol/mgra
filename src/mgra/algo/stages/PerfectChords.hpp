//
// Created by pavel on 1/19/16.
//

#ifndef MGRA_PERFECT_CHORDS_HPP
#define MGRA_PERFECT_CHORDS_HPP

namespace algo {

template<class graph_pack_t>
struct ProcessPerfectChords : public algo::AbsStage<graph_pack_t> {

    using mcolor_t = typename graph_pack_t::mcolor_t;

    using edge_t = typename graph_pack_t::edge_t;
    using arc_t = typename graph_pack_t::arc_t;
    using mularcs_t = typename graph_pack_t::mularcs_t;
    using twobreak_t = typename graph_pack_t::twobreak_t;
    using equiv_t = typename utility::equivalence<vertex_t>;


    explicit ProcessPerfectChords(size_t max_round)
            : algo::AbsStage<graph_pack_t>("perfect chords", "perfect_chords", round_stage_t, max_round) {
    }

    bool run(graph_pack_t &graph_pack) override;

private:
    std::map<vertex_t, size_t> determine_length_of_components(graph_pack_t &graph_pack,
                                                              std::map<vertex_t, std::set<vertex_t> > const &classes,
                                                              mcolor_t const &color);

    size_t get_length_of_component(graph_pack_t & graph_pack, vertex_t const & start, mcolor_t const &color);
private:
    DECL_LOGGER("ProcessPerfectChords");
};

template<class graph_pack_t>
bool ProcessPerfectChords<graph_pack_t>::run(graph_pack_t &graph_pack) {
    bool isChanged = false;
    size_t number_rear = 0; // number of rearrangements
    auto medians = graph_pack.multicolors.get_median_colors();
    assert(medians.size() == 1);

    mcolor_t left, right, parent;
    std::tie(left, right, parent) = medians[0];
    std::cerr << cfg::get().mcolor_to_name(left) << " " << cfg::get().mcolor_to_name(right) << " " <<
                 cfg::get().mcolor_to_name(parent) << std::endl;
    assert(left == right);

    do {
        number_rear = 0;

        equiv_t connected_components = graph_pack.split_on_components_with_color(parent);
        std::map<vertex_t, std::set<vertex_t> > const &classes = connected_components.get_eclasses<std::set<vertex_t>>();
        auto size_components = determine_length_of_components(graph_pack, classes, left);

        /*for (auto info : size_components)  {
            std::cerr << info.first << " l(e) = " << info.second << " l(v) = " << classes.find(info.first)->second.size() << std::endl;
        }*/

        for (auto const & info : classes) {
            if (size_components.find(info.first)->second % 2 == 0) {
                for (auto const & v : info.second) {
                    auto mularcs_v = graph_pack.get_all_adjacent_multiedges(v);

                    assert(mularcs_v.size() <= 3);

                    if (mularcs_v.size() != 3) {
                        continue;
                    }
                    //assert(mularcs_v.size() == 3);

                    vertex_t u = mularcs_v.get_vertex(parent);
                    assert(!u.empty());
                    if (u == Infty) {
                        continue;
                    }
                    auto const & mularcs_u = graph_pack.get_all_adjacent_multiedges(u);
                    mularcs_v.erase(u);
                    assert(mularcs_v.size() == 2);

                    if (size_components.find(connected_components[u])->second % 2 != 0) {
                        continue;
                    }

                    vertex_t const & p = mularcs_v.cbegin()->first;
                    if (p == Infty) {
                        continue;
                    }

                    if (graph_pack.get_all_adjacent_multiedges(p).size() != 3) {
                        continue;
                    }

                    vertex_t const & q = graph_pack.get_all_adjacent_multiedges(p).get_vertex(parent);
                    assert(!q.empty());

                    vertex_t const & w = mularcs_v.crbegin()->first;
                    if (w == Infty) {
                        continue;
                    }
                    if (graph_pack.get_all_adjacent_multiedges(w).size() != 3) {
                        continue;
                    }
                    vertex_t const & r = graph_pack.get_all_adjacent_multiedges(w).get_vertex(parent);
                    assert(!r.empty());

                    if (connected_components[u] == connected_components[q] &&
                        connected_components[u] != connected_components[r]) {
                        if (mularcs_u.defined(q)) {
                            ++number_rear;
                            twobreak_t twobreak(p, q, v, u, parent);
                            graph_pack.apply(twobreak);
                        }
                    } else if (connected_components[u] != connected_components[q] &&
                               connected_components[u] == connected_components[r]) {
                        if (mularcs_u.defined(r)) {
                            ++number_rear;
                            twobreak_t twobreak(p, q, w, r, parent);
                            graph_pack.apply(twobreak);
                        }
                    } else if (connected_components[u] == connected_components[q] &&
                               connected_components[u] == connected_components[r]) {
                        if (mularcs_u.defined(q) && !mularcs_u.defined(r)) {
                            ++number_rear;
                            twobreak_t twobreak(p, q, v, u, parent);
                            graph_pack.apply(twobreak);
                        } else if (!mularcs_u.defined(q) && mularcs_u.defined(r)) {
                            ++number_rear;
                            twobreak_t twobreak(p, q, w, r, parent);
                            graph_pack.apply(twobreak);
                        }
                    }

                    if (number_rear != 0) {
                        break;
                    }
                }

                if (number_rear != 0) {
                    break;
                }
            }
        }

        if (number_rear != 0) {
            std::cerr << "PROCESS HORD MODIFY SOMETHING" << std::endl;
            isChanged = true;
        }
    } while (number_rear > 0);

    return isChanged;
}

template<class graph_pack_t>
std::map<vertex_t, size_t> ProcessPerfectChords<graph_pack_t>::determine_length_of_components(graph_pack_t & graph_pack,
        std::map<vertex_t, std::set<vertex_t> > const &classes, mcolor_t const &color) {
    std::map<vertex_t, size_t> result;

    for (auto const & elem : classes) {
        size_t length = get_length_of_component(graph_pack, elem.first, color);
        result.insert(std::make_pair(elem.first, length));
    }

    return result;
}

template<class graph_pack_t>
size_t ProcessPerfectChords<graph_pack_t>::get_length_of_component(graph_pack_t & graph_pack, vertex_t const &start, mcolor_t const &color) {
    size_t length = 0;
    std::unordered_set<vertex_t> visited({start});
    std::stack<vertex_t> vertices_to_visit;


    for (auto const &edge: graph_pack.get_all_adjacent_multiedges(start)) {
        if (edge.second == color) {
            vertices_to_visit.push(edge.first);
        }
    }

    while (!vertices_to_visit.empty()) {
        auto const & v = vertices_to_visit.top();
        vertices_to_visit.pop();

        if (visited.count(v) != 0) {
            ++length;
            break;
        }

        ++length;
        visited.insert(v);

        auto const & mularcs = graph_pack.get_all_adjacent_multiedges(v);
        for (auto const &edge: mularcs) {
            if (edge.first != Infty && visited.count(edge.first) == 0) {
                if (edge.second == color) {
                    vertices_to_visit.push(edge.first);
                }
            }
        }
    }

    return length;
}

}

#endif //MGRA_PERFECT_CHORDS_HPP
