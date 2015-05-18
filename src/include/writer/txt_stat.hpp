#ifndef TXT_STATS_HPP
#define TXT_STATS_HPP

#include "io/path_helper.hpp"
#include "graph/graph_pack.hpp"
#include "structures/config_struct.hpp"
#include "reader.h"

namespace writer {

  template <class statistics_t>
  struct TXT_statistics {

    TXT_statistics() {
    }

    explicit TXT_statistics(std::string const& path_to_file)
        : ofstat(path_to_file) {
    }

    TXT_statistics(std::string const& path_to_file, std::string const& name_file)
        : ofstat(path::append_path(path_to_file, name_file)) {
    }

    void open(std::string const& path_to_file) {
      ofstat.open(path::append_path(path_to_file, "stats.tex"));
      print_opening_tex();
    }

    ~TXT_statistics() {
      if (ofstat.is_open()) {
        ofstat << "\\end{document}" << std::endl;
      }
    }

    void print_all_statistics(size_t stage, statistics_t const& info);

    void print_history_statistics(statistics_t const& info);

  private:
    void print_vertex_statistics(statistics_t const& info);

    void print_complete_edges(statistics_t const& info);

    void print_connected_components(statistics_t const& info);

    void print_rear_characters(statistics_t const& info);

    void print_indel_statistics(statistics_t const& info);

    void print_start_table(size_t count_column) {
      ofstat << "\\begin{table}[h]" << std::endl;
      ofstat << "\\centering \\begin{tabular}{|c|";
      for (size_t i = 0; i < count_column - 1; ++i) {
        ofstat << "c|";
      }
      ofstat << "}" << std::endl;
      ofstat << "\\hline" << std::endl;
    }

    void print_close_table() {
      ofstat << "\\hline" << std::endl;
      ofstat << "\\end{tabular}" << std::endl;
      ofstat << "\\end{table}" << std::endl;
      ofstat << std::endl;
    }

    void print_opening_tex() {
      ofstat << "\\documentclass[a4paper]{article}" << std::endl;
      ofstat << "\\usepackage[english]{babel}" << std::endl;
      ofstat << "\\usepackage[utf8]{inputenc}" << std::endl;
      ofstat << "\\usepackage{amsmath}" << std::endl;
      ofstat << "\\usepackage{graphicx}" << std::endl;
      ofstat << "\\title{MGRA output}" << std::endl;
      ofstat << "\\author{MGRA}" << std::endl;
      ofstat << "\\begin{document}" << std::endl;
    }

    void line_break() {
      ofstat << "\\\\" << std::endl;
    }

    std::ofstream ofstat;
  };

  template <class statistics_t>
  void TXT_statistics<statistics_t>::print_all_statistics(size_t stage, statistics_t const& info) {
    ofstat << "After " << stage << " stage: ";
    line_break();
    print_vertex_statistics(info);
    print_complete_edges(info);
    print_connected_components(info);
    print_rear_characters(info);
    print_indel_statistics(info);
  }

  template <class statistics_t>
  void TXT_statistics<statistics_t>::print_history_statistics(statistics_t const& info) {
    using mcolor_t = typename statistics_t::mcolor_type;
    auto print_lambda = [&](std::string const& name, std::map<mcolor_t, size_t> const& opers) -> void {
      ofstat << std::endl << name << opers.size() << " " << std::endl;
      for (auto const& event : opers) {
        ofstat << cfg::get().mcolor_to_name(event.first) << "\t" << event.second << std::endl;
      }
      ofstat << std::endl;
    };

    print_lambda("Total number of 2-breaks:", info.number_twobreaks);
    print_lambda("Total number of insertion events:", info.number_insertions);
    print_lambda("Total number of deletion events:", info.number_deletions);
  }

  template <class statistics_t>
  void TXT_statistics<statistics_t>::print_vertex_statistics(statistics_t const& info) {
    ofstat << "... Vertices: " << info.vertex_statistics[0];
    line_break();
    ofstat << "... Duplication vertex: " << info.vertex_statistics[1];
    line_break();
    ofstat << "... Insertion/deletion vertex: " << info.vertex_statistics[2];
    line_break();
    ofstat << "... Count self loop: " << info.vertex_statistics[3];
    line_break();
  }

  template <class statistics_t>
  void TXT_statistics<statistics_t>::print_complete_edges(statistics_t const& info) {
    ofstat << "... complete multiedges:";
    for (auto const& edge : info.complete_edges) {
      ofstat << " " << edge.first << "~" << edge.second;
    }
    ofstat << "\t(total: " << info.complete_edges.size() << ")";
    line_break();
  }

  template <class statistics_t>
  void TXT_statistics<statistics_t>::print_connected_components(statistics_t const& info) {
    ofstat << "... connected components:";
    for (auto const& sc : info.size_component_to_count) {
      ofstat << " " << "$" << sc.first << "^{\\wedge}" << sc.second << "$";
    }
    line_break();
  }

  template <class statistics_t>
  void TXT_statistics<statistics_t>::print_rear_characters(statistics_t const& info) {
    using mcolor_t = typename statistics_t::mcolor_type;

    ofstat << std::endl << "% Rearrangement characters:" << std::endl << std::endl;
    print_start_table(6);
    ofstat <<
    "Multicolors & multiedges & simple vertices & simple multiedges & simple paths+cycles & irreg. multiedges\\\\" <<
    std::endl;
    ofstat << "\\hline" << std::endl;

    auto calc_value = [](std::map<mcolor_t, size_t> const& where, mcolor_t const& what) -> size_t {
      if (where.find(what) != where.end()) {
        return where.find(what)->second;
      }
      return 0;
    };

    std::multimap<size_t, std::string> answer;
    for (auto const& edge : info.multiedges_count) {
      mcolor_t const& current = edge.first.second;
      size_t res = edge.second / 2;
      size_t paths = calc_value(info.simple_vertices_count, edge.first.second) - \
          (
          calc_value(info.simple_multiedges_count, edge.first.second) + \
          calc_value(info.simple_multiedges_count, edge.first.first)\
) - \
          calc_value(info.simple_vertices_alone_count, edge.first.second) - \
          calc_value(info.special_cycle_count, edge.first.second);
      size_t cycles = calc_value(info.simple_cycle_count, current) + \
                    calc_value(info.special_cycle_count, current);


      std::ostringstream os;

      os << "{" << cfg::get().mcolor_to_name(edge.first.second) << " + "
      << cfg::get().mcolor_to_name(edge.first.first) << "} & "
      // multiedges
      << res << " & "
      // simple vertices
      << calc_value(info.simple_vertices_count, current) << " & "
      // simple multiedges
      << calc_value(info.simple_multiedges_count, current) << " + "
      << calc_value(info.simple_multiedges_count, edge.first.first)
      << " = " <<
      calc_value(info.simple_multiedges_count, current) + calc_value(info.simple_multiedges_count, edge.first.first) <<
      " & "
      // simple paths + cycles
      << paths << " + " << cycles << " = " << paths + cycles << " & "
      // irregular multiedges
      << calc_value(info.irrer_multiedges_count, current) << " + "
      << calc_value(info.irrer_multiedges_count, edge.first.first) << " = "
      << calc_value(info.irrer_multiedges_count, current) + calc_value(info.irrer_multiedges_count, edge.first.first);

      answer.insert(std::make_pair(res, os.str()));
    }

    for (auto im = answer.rbegin(); im != answer.rend(); ++im) {
      ofstat << im->second << "\\\\" << std::endl;
    }

    print_close_table();
  }

  template <class statistics_t>
  void TXT_statistics<statistics_t>::print_indel_statistics(statistics_t const& info) {
    ofstat << std::endl << "% Insertion/Deletion characters: " << std::endl << std::endl;
    print_start_table(2);

    ofstat << "insert multicolor $Q + \\bar{Q}$ & number edges \\\\" << std::endl;
    ofstat << "\\hline" << std::endl;

    for (auto const& indel : info.complement_indel_stats) {
      ofstat << cfg::get().mcolor_to_name(indel.first.first) << " + "
      << cfg::get().mcolor_to_name(indel.first.second) << " & " << indel.second << "\\\\" << std::endl;
    }

    print_close_table();
  }
}

#endif
