#ifndef AS4_STAGE_HPP
#define AS4_STAGE_HPP

namespace algo {

template<class graph_pack_t>
struct ProcessAS4Patterns : public algo::AbsStage<graph_pack_t> {

    using mcolor_t = typename graph_pack_t::mcolor_type;
    using edge_t = typename graph_pack_t::edge_t;
    using arc_t = typename graph_pack_t::arc_t;
    using mularcs_t = typename graph_pack_t::mularcs_t;
    using twobreak_t = typename graph_pack_t::twobreak_t;


    explicit ProcessAS4Patterns(size_t max_round = 1)
        : AbsStage<graph_pack_t>("Process AS4 patterns", "as4_patterns", round_stage_t, max_round)
    {
    }

    bool run(graph_pack_t& graph_pack) override;

private:
    boost::optional<vertex_t> get_connection(vertex_t& v1, vertex_t& v2,
            graph_pack_t& graph_pack, std::array<mcolor_t, 3>& colors);
    bool is_connected(vertex_t& v1, vertex_t& v2,
            graph_pack_t& graph_pack, std::array<mcolor_t, 3>& colors);
    bool process_as4_patterns(graph_pack_t& graph_pack);
};

template <class graph_pack_t>
bool ProcessAS4Patterns<graph_pack_t>::is_connected(vertex_t &v1, vertex_t &v2,
        graph_pack_t &graph_pack, std::array<mcolor_t, 3> &colors)
{
    for (size_t idx = 0; idx < 3; ++idx)
    {
        vertex_t cand = graph_pack.get_all_adjacent_multiedges(v1).get_vertex(colors[idx]);
        if (cand == v2) return true;
    }
    return false;
}
template <class graph_pack_t>
boost::optional<vertex_t> ProcessAS4Patterns<graph_pack_t>::get_connection(vertex_t& v1, vertex_t& v2,
        graph_pack_t& graph_pack, std::array<mcolor_t, 3>& colors)
{
    boost::optional<vertex_t> res;
    for (size_t idx = 0; idx < 3; ++idx)
    {
        vertex_t candidate = graph_pack.get_all_adjacent_multiedges(v1).get_vertex(colors[idx]);
        if (candidate == Infty || graph_pack.graph.degree_vertex(candidate) != 3) continue;
        for (size_t j = 1; j < 3; ++j)
        {
            vertex_t match = graph_pack.get_all_adjacent_multiedges(candidate).get_vertex(colors[(idx+j) % 3]);
            if (match == v2)
            {
                res = candidate;
                return res;
            }
        }
    }
    return res;
}
template<class graph_pack_t>
bool ProcessAS4Patterns<graph_pack_t>::process_as4_patterns(graph_pack_t& graph_pack) {
    TRACE("Process AS4 patterns ")
    bool isChanged = false;
    auto medians = graph_pack.multicolors.get_medians_colors();
    std::vector<edge_t> ancestor_adjacency;
    std::unordered_set<vertex_t> in_ancestor;

    for(auto const & median_color : medians) {
        mcolor_t left, right, parent;
        std::tie(left, right, parent) = median_color;
        std::string sleft = cfg::get().mcolor_to_name(left);
        std::string sright = cfg::get().mcolor_to_name(right);
        std::string sparent = cfg::get().mcolor_to_name(parent);
        TRACE("Start to work with median " << sleft << " " << sright << " " << sparent)

        std::array<mcolor_t, 3> colors;
        colors[0] = left; colors[1] = right; colors[2] = parent;

        for (size_t ind_color = 0; ind_color < 3; ++ind_color) {
            mcolor_t c0 = colors[ind_color]; mcolor_t c1 = colors[(ind_color + 1) % 3]; mcolor_t c2 = colors[(ind_color + 2) % 3];
            mcolor_t c01(c0, c1, mcolor_t::Union);
            std::unordered_set<vertex_t> marked = in_ancestor;

            for (vertex_t const & i : graph_pack.graph) {
                if (graph_pack.graph.degree_vertex(i) != 3) continue;

                vertex_t j = graph_pack.get_all_adjacent_multiedges(i).get_vertex(c0);
                if (!j.empty())
                {
                    if (j == Infty || graph_pack.graph.degree_vertex(j) != 3) continue;
                    if (marked.count(i) != 0 || marked.count(j) != 0) continue;
                    
                    vertex_t p_u_l = graph_pack.get_all_adjacent_multiedges(i).get_vertex(c1);
                    vertex_t p_d_l = graph_pack.get_all_adjacent_multiedges(i).get_vertex(c2);
                    vertex_t p_u_r = graph_pack.get_all_adjacent_multiedges(j).get_vertex(c2);
                    vertex_t p_d_r = graph_pack.get_all_adjacent_multiedges(j).get_vertex(c1);

                    TRACE("See on " << i << " " << j << " " << p_u_l << " " << p_d_l << " " << p_u_r << " " << p_d_r)
                    auto reg_check = [&] (vertex_t v) { return !v.empty() && marked.count(v) == 0 &&
                            v != Infty && graph_pack.graph.degree_vertex(v) == 3; };
                    if (reg_check(p_d_l) && reg_check(p_d_r) && reg_check(p_u_r) && reg_check(p_u_l)) {
                        vertex_t m0, m1, m2;
                        if (p_u_l == p_u_r)
                        {
                            vertex_t v1 = graph_pack.get_all_adjacent_multiedges(p_d_l).get_vertex(c1);
                            vertex_t v2 = graph_pack.get_all_adjacent_multiedges(p_d_r).get_vertex(c2);
                            if (v1 == v2)
                                //upper envelope case
                            {
                                m0 = v1;
                                m1 = graph_pack.get_all_adjacent_multiedges(p_d_l).get_vertex(c0);
                                m2 = graph_pack.get_all_adjacent_multiedges(p_d_r).get_vertex(c0);
                                std::unordered_set<vertex_t> counter;
                                counter.insert(i); counter.insert(p_d_l); counter.insert(j); counter.insert(p_d_r);
                                counter.insert(p_u_l); counter.insert(m0); counter.insert(m1); counter.insert(m2);
                                if (counter.size() != 8) continue;
                                if (!reg_check(m0) || !reg_check(m1) || !reg_check(m2)) continue;
                                if (graph_pack.get_all_adjacent_multiedges(m1).get_vertex(c1) == m2
                                        || graph_pack.get_all_adjacent_multiedges(m1).get_vertex(c2) == m2)
                                {
                                    isChanged = true;
                                    TRACE("TYPE 1-UP")
                                    ancestor_adjacency.push_back(edge_t(i, p_d_l));
                                    ancestor_adjacency.push_back(edge_t(j, p_d_r));
                                    ancestor_adjacency.push_back(edge_t(p_u_l, m0));
                                    ancestor_adjacency.push_back(edge_t(m1, m2));
                                    in_ancestor.insert(i); in_ancestor.insert(j); in_ancestor.insert(p_d_l); in_ancestor.insert(p_d_r);
                                    marked.insert(i); marked.insert(j); marked.insert(p_d_l); marked.insert(p_d_r);
                                    in_ancestor.insert(p_u_l); in_ancestor.insert(m0); in_ancestor.insert(m1); in_ancestor.insert(m2);
                                    marked.insert(p_u_l); marked.insert(m0); marked.insert(m1); marked.insert(m2);
                                }
                            }
                        }
                        else if (p_d_l == p_d_r)
                        {
                            vertex_t v1 = graph_pack.get_all_adjacent_multiedges(p_u_l).get_vertex(c2);
                            vertex_t v2 = graph_pack.get_all_adjacent_multiedges(p_u_r).get_vertex(c1);
                            if (v1 != v2) continue;
                            m0 = v1;
                            m1 = graph_pack.get_all_adjacent_multiedges(p_u_l).get_vertex(c0);
                            m2 = graph_pack.get_all_adjacent_multiedges(p_u_r).get_vertex(c0);
                            std::unordered_set<vertex_t> counter;
                            counter.insert(i); counter.insert(p_u_l); counter.insert(j); counter.insert(p_u_r);
                            counter.insert(p_d_l); counter.insert(m0); counter.insert(m1); counter.insert(m2);
                            if (counter.size() != 8) continue;
                            if (!reg_check(m0) || !reg_check(m1) || !reg_check(m2)) continue;
                            if (graph_pack.get_all_adjacent_multiedges(m1).get_vertex(c1) == m2
                                    || graph_pack.get_all_adjacent_multiedges(m1).get_vertex(c2) == m2)
                            {
                                isChanged = true;
                                TRACE("TYPE 1-DOWN")
                                ancestor_adjacency.push_back(edge_t(i, p_u_l));
                                ancestor_adjacency.push_back(edge_t(j, p_u_r));
                                ancestor_adjacency.push_back(edge_t(p_d_l, m0));
                                ancestor_adjacency.push_back(edge_t(m1, m2));
                                in_ancestor.insert(i); in_ancestor.insert(j); in_ancestor.insert(p_u_l); in_ancestor.insert(p_u_r);
                                marked.insert(i); marked.insert(j); marked.insert(p_u_l); marked.insert(p_u_r);
                                in_ancestor.insert(p_d_l); in_ancestor.insert(m0); in_ancestor.insert(m1); in_ancestor.insert(m2);
                                marked.insert(p_d_l); marked.insert(m0); marked.insert(m1); marked.insert(m2);
                            }
                        }
                        else if (graph_pack.get_all_adjacent_multiedges(p_u_l).get_vertex(c0) == p_d_l &&
                                graph_pack.get_all_adjacent_multiedges(p_u_r).get_vertex(c0) == p_d_r)
                        {
                            vertex_t n1 = graph_pack.get_all_adjacent_multiedges(p_u_l).get_vertex(c2);
                            vertex_t n2 = graph_pack.get_all_adjacent_multiedges(p_d_l).get_vertex(c1);
                            vertex_t n3 = graph_pack.get_all_adjacent_multiedges(p_d_r).get_vertex(c2);
                            vertex_t n4 = graph_pack.get_all_adjacent_multiedges(p_u_r).get_vertex(c1);
                            vertex_t e1, e2;
                            if (is_connected(n2, n4, graph_pack, colors))
                            {
                                e1 = n2;
                                e2 = n4;
                            }
                            else if (is_connected(n1, n3, graph_pack, colors))
                            {
                                e1 = n1;
                                e2 = n3;
                            }
                            else
                            {
                                continue;
                            }
                            std::unordered_set<vertex_t> counter;
                            counter.insert(i); counter.insert(j); counter.insert(e1); counter.insert(e2);
                            counter.insert(p_d_l); counter.insert(p_u_r); counter.insert(p_u_l); counter.insert(p_d_r);
                            if (counter.size() != 8) continue;
                            else
                            {
                                TRACE("EARS-CASE")
                                ancestor_adjacency.push_back(edge_t(i, j));
                                ancestor_adjacency.push_back(edge_t(e1, e2));
                                ancestor_adjacency.push_back(edge_t(p_d_l, p_u_r));
                                ancestor_adjacency.push_back(edge_t(p_d_r, p_u_l));
                                in_ancestor.insert(i); in_ancestor.insert(j); in_ancestor.insert(p_d_l); in_ancestor.insert(p_d_r);
                                marked.insert(i); marked.insert(j); marked.insert(p_d_l); marked.insert(p_d_r);
                                in_ancestor.insert(p_u_l); in_ancestor.insert(p_u_r); in_ancestor.insert(e1); in_ancestor.insert(e2);
                                marked.insert(p_u_l); marked.insert(p_u_r); marked.insert(e1); marked.insert(e2);
                            }
                        }
                        else
                        {
                            // second type branch
                            auto up = get_connection(p_u_l, p_u_r, graph_pack, colors);
                            auto down = get_connection(p_d_l, p_d_r, graph_pack, colors);
                            if (up && down)
                            {

                                vertex_t p_u_m = *up;
                                vertex_t p_d_m = *down;
                                if (!reg_check(p_u_m) || !reg_check(p_d_m)) continue;
                                std::unordered_set<vertex_t> counter;
                                counter.insert(i); counter.insert(p_d_l); counter.insert(j); counter.insert(p_u_l);
                                counter.insert(p_d_r); counter.insert(p_d_m); counter.insert(p_u_r); counter.insert(p_u_m);
                                if (counter.size() != 8) continue;
                                vertex_t v1 = graph_pack.get_all_adjacent_multiedges(p_u_l).get_vertex(c0);
                                if (!v1.empty() && v1 == p_d_l)
                                {
                                    TRACE("TYPE 2-LEFT")
                                    isChanged = true;
                                    ancestor_adjacency.push_back(edge_t(i, p_d_l));
                                    ancestor_adjacency.push_back(edge_t(p_d_r, p_d_m));
                                    ancestor_adjacency.push_back(edge_t(j, p_u_l));
                                    ancestor_adjacency.push_back(edge_t(p_u_r, p_u_m));
                                    in_ancestor.insert(i); in_ancestor.insert(p_d_l); in_ancestor.insert(j); in_ancestor.insert(p_u_l);
                                    marked.insert(i); marked.insert(p_d_l); marked.insert(j); marked.insert(p_u_l);
                                    in_ancestor.insert(p_d_r); in_ancestor.insert(p_d_m); in_ancestor.insert(p_u_r); in_ancestor.insert(p_u_m);
                                    marked.insert(p_d_r); marked.insert(p_d_m); marked.insert(p_u_r); marked.insert(p_u_m);
                                    continue;
                                }
                                vertex_t v2 = graph_pack.get_all_adjacent_multiedges(p_u_r).get_vertex(c0);
                                if (!v2.empty() && v2 == p_d_r)
                                {
                                    TRACE("TYPE 2-RIGHT")
                                    isChanged = true;
                                    ancestor_adjacency.push_back(edge_t(j, p_d_r));
                                    ancestor_adjacency.push_back(edge_t(p_d_l, p_d_m));
                                    ancestor_adjacency.push_back(edge_t(i, p_u_r));
                                    ancestor_adjacency.push_back(edge_t(p_u_l, p_u_m));
                                    in_ancestor.insert(i); in_ancestor.insert(p_d_l); in_ancestor.insert(j); in_ancestor.insert(p_u_l);
                                    marked.insert(i); marked.insert(p_d_l); marked.insert(j); marked.insert(p_u_l);
                                    in_ancestor.insert(p_d_r); in_ancestor.insert(p_d_m); in_ancestor.insert(p_u_r); in_ancestor.insert(p_u_m);
                                    marked.insert(p_d_r); marked.insert(p_d_m); marked.insert(p_u_r); marked.insert(p_u_m);
                                    continue;
                                }
                                vertex_t vc0 = graph_pack.get_all_adjacent_multiedges(p_u_m).get_vertex(c0);
                                vertex_t vc1 = graph_pack.get_all_adjacent_multiedges(p_u_m).get_vertex(c1);
                                vertex_t vc2 = graph_pack.get_all_adjacent_multiedges(p_u_m).get_vertex(c2);
                                if (vc0 == p_d_m || vc1 == p_d_m || vc2 == p_d_m)
                                {
                                    TRACE("TYPE 2-CENTER")
                                    isChanged = true;
                                    ancestor_adjacency.push_back(edge_t(i, j));
                                    ancestor_adjacency.push_back(edge_t(p_u_m, p_d_m));
                                    ancestor_adjacency.push_back(edge_t(p_u_l, p_d_r));
                                    ancestor_adjacency.push_back(edge_t(p_u_r, p_d_l));
                                    in_ancestor.insert(i); in_ancestor.insert(p_d_l); in_ancestor.insert(j); in_ancestor.insert(p_u_l);
                                    marked.insert(i); marked.insert(p_d_l); marked.insert(j); marked.insert(p_u_l);
                                    in_ancestor.insert(p_d_r); in_ancestor.insert(p_d_m); in_ancestor.insert(p_u_r); in_ancestor.insert(p_u_m);
                                    marked.insert(p_d_r); marked.insert(p_d_m); marked.insert(p_u_r); marked.insert(p_u_m);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (edge_t const & edge : ancestor_adjacency) {
        if (graph_pack.get_all_multicolor_edge(edge.first, edge.second) != graph_pack.multicolors.get_complete_color()) {
            std::string central_color = cfg::get().mcolor_to_name(graph_pack.get_all_multicolor_edge(edge.first, edge.second));
            TRACE("Apply twobreak on " << edge.first << " " << edge.second << " " << central_color)

            mularcs_t mularcs_first = graph_pack.get_all_adjacent_multiedges(edge.first);
            mularcs_first.erase(edge.second);
            mularcs_t mularcs_second = graph_pack.get_all_adjacent_multiedges(edge.second);
            mularcs_second.erase(edge.first);

            for (arc_t const & arc : mularcs_first) {
                if (graph_pack.multicolors.is_vec_T_consistent_color(arc.second)) {
                    isChanged = true;
                    vertex_t v = mularcs_second.get_vertex(arc.second);
                    assert(!v.empty());
                    graph_pack.apply(twobreak_t(edge.first, arc.first, edge.second, v, arc.second));
                }
            }
        }
    }
    return isChanged;
}

template<class graph_pack_t>
bool ProcessAS4Patterns<graph_pack_t>::run(graph_pack_t& graph_pack) {
    bool isChanged = false;
    bool temp = true;

    while(temp) {
        temp = false;
        temp = process_as4_patterns(graph_pack);
        if (temp) {
            isChanged = true;
        }
    }
    TRACE("Finish process special vec{T}-consistent subgraphs " << isChanged)
    return isChanged;
}

}
#endif //_MGRA_AS4_HPP_