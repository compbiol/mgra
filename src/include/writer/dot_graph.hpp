#ifndef GRAPHDOTS_HPP
#define GRAPHDOTS_HPP

namespace writer {

  template <class graph_pack_t>
  struct GraphDot {
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using edge_t = typename graph_pack_t::edge_t;
    using mularcs_t = typename graph_pack_t::mularcs_t;
    using tree_t = typename structure::BinaryTree<mcolor_t>;
    using node_t = typename tree_t::colored_node_t;

    GraphDot()
        : infinity_verteces(0) {
    }

    void open(std::string const& path, std::string const& graphname = "stage") {
      m_path = path;
      m_graphname = graphname;
    }

    /**
    * Save multiply breakpoint graph in .dot file
    */
    void save_bp_graph(graph_pack_t const& graph, size_t stage);

    /**
    * Save subtrees in .dot file with legend colors.
    */
    void save_subtrees();

    void save_subtrees(std::vector <tree_t> const& trees);

    void save_subtrees(std::ostream& output,
        std::vector <tree_t> const& trees);

  private:
    std::string get_branches_from_tree(std::shared_ptr<const node_t> const& current, std::vector <std::string>& info) const;

    // Save multiedge in dot file
    void draw_multiedge(std::ofstream& dot, vertex_t const& x, mcolor_t const& color, vertex_t const& y,
        std::string const& proporties, bool is_label_proporties = false);

    std::string m_path;
    std::string m_graphname;
    int infinity_verteces;
  };
}

#include "dot_graph_impl.hpp"

#endif

