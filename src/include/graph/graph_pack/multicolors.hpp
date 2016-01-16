#ifndef MULTICOLORS_HPP
#define MULTICOLORS_HPP

namespace graph {

namespace graph_pack {

//TODO introduce T-solid color.

template<class mcolor_type>
struct Multicolors {
    using mcolor_t = mcolor_type;
    using mcolor_ptr = typename mcolor_t::mcolor_ptr;
    using citer = typename std::set<mcolor_t>::const_iterator;

    using genome_number_t = std::unordered_map<std::string, size_t>;
    using colors_median_t = std::tuple<mcolor_t, mcolor_t, mcolor_t>;
    using vector_colors_median_t = std::vector<colors_median_t>;

    template<class tree_t>
    Multicolors(tree_t const & tree, genome_number_t const & genome_number, std::string const & target_name);

    /*
     * Todo
     */
    mcolor_t const & get_complement_color(mcolor_t const &color) {
        if (complement_colors.count(color) == 0) {
            assert(complete_color.includes(color));
            mcolor_t temp(complete_color, color, mcolor_t::Difference);
            complement_colors.insert(std::make_pair(color, temp));
            complement_colors.insert(std::make_pair(temp, color));
        }
        return complement_colors.find(color)->second;
    }

    /*
     * Split input multiolor COLOR on T-consistent multicolors according number of splits.
     * If multicolor split to more than NUMBER_SPLITS value return input color.
     * If less returned set of T-consistent color.
     */
    std::set<mcolor_t> split_color_on_tc_color(mcolor_t const &color, size_t number_splits = 1) const{
        if (is_T_consistent_color(color) || (number_splits == 1)) {
            return std::set<mcolor_t>({color});
        } else {
            std::set<mcolor_t> answer = split_color(color, T_consistent_colors);
            if (answer.size() > number_splits) {
                return std::set<mcolor_t>({color});
            }
            return answer;
        }
    }

    /*
     * Split input multiolor COLOR on vec{T}-consistent multicolors.
     */
    std::set<mcolor_t> split_color_on_vtc_color(mcolor_t const &color) const {
        if (is_vec_T_consistent_color(color)) {
            return std::set<mcolor_t>({color});
        } else {
            return split_color(color, vec_T_consistent_colors);
        }
    }

    /*
     * TODO
     */
    mcolor_t get_min_addit_color_for_tc(mcolor_t const &color) {
        assert(complete_color.includes(color));

        mcolor_t min_color = get_complete_color();
        for (auto col = cbegin_T_consistent_color(); col != cend_T_consistent_color(); ++col) {
            if (col->includes(color)) {
                mcolor_t diff_color(*col, color, mcolor_t::Difference);
                if (diff_color.size() < min_color.size()) {
                    min_color = diff_color;
                }
            }
        }
        return min_color;
    }

    bool is_T_consistent_color(mcolor_t const &color) const {
        return (T_consistent_colors.count(color) > 0);
    }

    bool is_vec_T_consistent_color(mcolor_t const &color) const {
        return (vec_T_consistent_colors.count(color) > 0);
    }

    bool is_T_solid_color(mcolor_t const &color) const {
        return (T_solid_colors.count(color) > 0);
    }

    DECLARE_GETTER(vector_colors_median_t, median_colors, median_colors)

    DECLARE_GETTER(mcolor_t const &, complete_color, complete_color)

    DECLARE_GETTER(mcolor_t const &, remove_color, root_color)

    DECLARE_CONST_ITERATOR(citer, T_consistent_colors, cbegin_T_consistent_color, cbegin)

    DECLARE_CONST_ITERATOR(citer, T_consistent_colors, cend_T_consistent_color, cend)

    DECLARE_CONST_ITERATOR(citer, T_solid_colors, cbegin_T_solid_color, cbegin)

    DECLARE_CONST_ITERATOR(citer, T_solid_colors, cend_T_solid_color, cend)

    DECLARE_CONST_ITERATOR(citer, vec_T_consistent_colors, cbegin_vec_T_consistent_color, cbegin)

    DECLARE_CONST_ITERATOR(citer, vec_T_consistent_colors, cend_vec_T_consistent_color, cend)

private:
    struct McolorCompBySize {
        bool operator()(mcolor_t const & first, mcolor_t const & second) const {
            return (first.size() > second.size());
        }
    };

    /*
     *
     */
    std::set<mcolor_t> split_color(mcolor_t const & input_color, std::set<mcolor_t> const & split_colors) const;

protected:
    mcolor_t complete_color;
    mcolor_t remove_color;

    std::set<mcolor_t> T_consistent_colors; // vec_T_consistent_colors U T_solid colors
    std::set<mcolor_t> T_solid_colors; //Colors are complement to corresponding subtree based on choose root
    std::set<mcolor_t> vec_T_consistent_colors; // Colors corresponding subtree based on choose root

    std::map<mcolor_t, mcolor_t> complement_colors;

    vector_colors_median_t median_colors;

    std::map<mcolor_t, std::string> multicolors_to_name;
    std::unordered_map<std::string, mcolor_t> name_to_multicolors;
};

template<class mcolor_t>
template<class phyltree_t>
Multicolors<mcolor_t>::Multicolors(phyltree_t const &phylotree, genome_number_t const &genome_number,
                                   std::string const &target_name) {
    assert(!genome_number.empty());

    // Get all vec{T}-consistent colors corresponding input [sub]tree[s].
    structure::phyl_tree::get_T_multicolors(phylotree, genome_number, multicolors_to_name, complete_color);
    for (auto const & info : multicolors_to_name) {
        if (info.first != complete_color) {
            vec_T_consistent_colors.insert(info.first);
            T_solid_colors.insert(get_complement_color(info.first));
        }
        assert(name_to_multicolors.count(info.second) == 0);
        name_to_multicolors.insert(std::make_pair(info.second, info.first));
    }

    // Get all median colors corresponding input [sub]tree[s].
    structure::phyl_tree::get_median_colors(phylotree, complete_color, name_to_multicolors, median_colors);

    //If target is empty we put root in nearest node.
    if (phylotree.get_root()->get_type() == phyltree_t::node_t::classic) {
        assert(name_to_multicolors.count(phylotree.get_root()->get_most_left_child()->get_data()) != 0);
        assert(name_to_multicolors.count(phylotree.get_root()->get_most_right_child()->get_data()) != 0);
        mcolor_t const &left_color = name_to_multicolors.find(
                phylotree.get_root()->get_most_left_child()->get_data())->second;
        mcolor_t const &right_color = name_to_multicolors.find(
                phylotree.get_root()->get_most_right_child()->get_data())->second;
        remove_color = (left_color.size() > right_color.size()) ? left_color : right_color;
    }
    vec_T_consistent_colors.erase(remove_color);
    T_solid_colors.erase(get_complement_color(remove_color));

    if (!target_name.empty()) {
        if (name_to_multicolors.find(target_name) != name_to_multicolors.cend()) {
            assert(target_name != phylotree.get_root()->get_data());
            ;
            //remove_color = name_to_multicolors.find(target_name)->second;
            //TODO think about target color
        } else {
            FATAL_ERROR("ERROR: name of target must have name equal to name of corresponding node")
        }
    }

    //init T consistent multicolors (vec{T}-consistent multicolors and addition).
    T_consistent_colors.insert(vec_T_consistent_colors.begin(), vec_T_consistent_colors.end());
    T_consistent_colors.insert(T_solid_colors.begin(), T_solid_colors.end());
}

template<class mcolor_t>
std::set<mcolor_t> Multicolors<mcolor_t>::split_color(mcolor_t const & input_color, std::set<mcolor_t> const & split_colors) const {
    assert(complete_color.includes(input_color));

    std::multiset<mcolor_t, McolorCompBySize> inter_colors;
    for (auto const & split_color: split_colors) {
        if (input_color.includes(split_color)) {
            inter_colors.insert(split_color);
        }
    }

    std::set<mcolor_t> answer;
    mcolor_t color = input_color;
    for (auto const & inter_color: inter_colors) {
        if (color.includes(inter_color)) {
            color = mcolor_t(color, inter_color, mcolor_t::Difference);
            answer.insert(inter_color);
        }
    }

    if (color.empty()) {
        return answer;
    } else {
        return std::set<mcolor_t>({input_color});
    }
}

}

/*
 * Split input vec{T}-consistent multiolor COLOR on vec{T}-consistent multicolors.
 */
/*std::set<mcolor_t> split_color_on_next_vtc_color(mcolor_t const &color) const {
    assert(is_vec_T_consistent_color(color));
    std::set<mcolor_t> answer = split_color(color, vec_T_consistent_colors);
    if (answer.size() == 1 && *answer.begin() == color) {
        answer.clear();
    }
    return answer;
}*/

/*
methods for split color. works only on non-wgd colors
utility::equivalence<size_t> equiv;
std::for_each(color.cbegin(), color.cend(), [&](std::pair<size_t, size_t> const &col) -> void {
    equiv.addrel(col.first, col.first);
});

for (auto const & split_color: split_colors) {
    mcolor_t inter_color(split_color, color, mcolor_t::Intersection);
    if (inter_color.size() >= 2 && inter_color.size() == split_color.size()) {
        std::for_each(inter_color.cbegin(), inter_color.cend(), [&](std::pair<size_t, size_t> const &col) -> void {
            equiv.addrel(col.first, inter_color.cbegin()->first);
        });
    }
}

equiv.update();
std::map<size_t, mcolor_t> const &classes = equiv.get_eclasses<mcolor_t>();
for (auto const &col : classes) {
    answer.insert(col.second);
}

return answer;*/

}
#endif

