#ifndef WOUTPUT_HPP_
#define WOUTPUT_HPP_ 

#include "defined.h"
#include "graph/estimate.h"
#include "writer/Wstats.h"
#include "writer/Wdots.h"


namespace writer { 

struct Woutput { 
	Woutput() {
	}

	Woutput(ProblemInstance<structure::Mcolor> const & cfg, fs::path const & out_dir
		, std::string const & stats_file, std::string const & graph_file, std::string const & colorsheme) 
	: m_cfg(cfg)
  	, m_out_dir(out_dir) 
	, m_debug_dir(out_dir / "debug")
	, m_genomes_dir(out_dir / "genomes")
	, m_transformation_dir(out_dir / "transformations")
	, m_out_stats(out_dir / stats_file)
	, m_out_dots(m_debug_dir, colorsheme, graph_file)
	{ 
	}

	inline void write_legend_dot() { 
		m_out_dots.write_legend_dot(m_cfg);
	}

	inline void write_info_after_stage(size_t stage, Statistics<mbgraph_with_history<structure::Mcolor> >& info
		, mbgraph_with_history<structure::Mcolor> const & graph) {
		m_out_stats.print_all_statistics(stage, info, graph);
		m_out_dots.save_dot(graph, m_cfg, stage);
	}

	inline void end_write_debug(mbgraph_with_history<structure::Mcolor> const & graph, edges_t const & bad_edges) {
		m_out_dots.save_dot(graph, m_cfg, 99);
		m_out_stats.print_history_statistics(graph, bad_edges);
	}

	bool initilizate() {
		auto creater_lambda = [](fs::path const & directory) -> bool {
			sys::error_code error; 
			if ( fs::exists(directory, error) ) {
				if ( !fs::is_directory(directory, error) ) {
					return false;
				}
			} else {
				bool flag = fs::create_directory(directory, error); 

				if ( !flag || error != 0 ) { 
					return false;
				}
			}
			return true; 
		};

		return creater_lambda(m_debug_dir) && creater_lambda(m_genomes_dir) && creater_lambda(m_transformation_dir);
	}


private: 
	ProblemInstance<structure::Mcolor> m_cfg; 
	fs::path m_out_dir;
	fs::path m_debug_dir;
	fs::path m_genomes_dir; 
	fs::path m_transformation_dir; 

	writer::Wstats m_out_stats; 
	writer::Wdots<mbgraph_with_history<structure::Mcolor>, ProblemInstance<structure::Mcolor> > m_out_dots;
};

}

#endif