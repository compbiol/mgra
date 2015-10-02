/* 
** Module: MGRA main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
*/

#include "command_line_parsing.hpp"

#include "algo/Algorithms.hpp"
#include "algo/recover_tree/recover_tree_task.hpp"

#include "writer/txt_genome.hpp"
#include "writer/txt_transform.hpp"

int main(int argc, char** argv) {
  if (parse_config_from_command_line(argc, argv)) {
    std::cerr << "ERROR: error while parsing command line arguments" << std::endl;
    return 1;
  }

  if (validate_application_config()) {
    std::cerr << "ERROR: error while validating config" << std::endl;
    return 1;
  }

  if (!organize_output_directory()) {
    std::cerr << "ERROR: problem with organize output directory "
        << cfg::get().out_path_directory << std::endl;
    return 1;
  }

  create_logger_from_config();

  //Reading problem configuration and genomes
  using genome_t = structure::Genome;
  using mcolor_t = structure::Mcolor;
  using graph_pack_t = GraphPack<mcolor_t>;

  INFO("Parse genomes file")
  std::vector<genome_t> genomes;
  if (cfg::get().block_file_type == infercars) {
    genomes = reader::read_infercars(cfg::get().blocks_file_path);
  } else if (cfg::get().block_file_type == grimm) {
    genomes = reader::read_grimm(cfg::get().blocks_file_path);
  }

  for (size_t i = 0; i < genomes.size(); ++i) {
    std::ostringstream out;
    out << "Download genome " << cfg::get().get_priority_name(i) << " with " << genomes[i].size() << " blocks.";
    INFO(out.str())
  }

  //Do job
  INFO("Start build graph")
  graph_pack_t graph_pack(genomes);
  INFO("End build graph")

  if (cfg::get().is_recover_tree) {
    algo::recover_tree_task(graph_pack);
  } else {
    {
      std::ostringstream out;
      out << "Determine " << graph_pack.multicolors.count_vec_T_consitent_color() << " \\vec{T}-consistent colors in tree:\n";
      for (auto id = graph_pack.multicolors.cbegin_vec_T_consistent_color(); id != graph_pack.multicolors.cend_vec_T_consistent_color(); ++id) {
        out << cfg::get().mcolor_to_name(*id) << " ";
      }
      INFO(out.str())
    }

    boost::optional<algo::RecoveredInformation<graph_pack_t>::AncestorInformation> result; 
    if (cfg::get().how_build == default_algo) { 
      result = algo::main_algorithm(graph_pack);
    } else if (cfg::get().how_build == wgd_algo) { 
      result = algo::wgd_algorithm(graph_pack);
    } 

    /*Save different output information in files*/
    if (result) { 
      using ancestor_information_t = algo::RecoveredInformation<graph_pack_t>::AncestorInformation;
      ancestor_information_t info = *result; 

      if (!info.transformations.empty()) {
        INFO("Save transformations in files.")

        writer::TXT_transformation<graph_pack_t> writer_transform(cfg::get().out_path_directory, graph_pack); 

        for (auto const & transformation : info.transformations) { 
          writer_transform.save_transformation(transformation.first, transformation.second);
          writer_transform.save_reverse_transformation(transformation.first, transformation.second);
        }    
      } 

      if (!info.genomes.empty()) {
        INFO("Save ancestor genomes in files.")
        writer::TXT_genome<genome_t> writer_genome(path::append_path(cfg::get().out_path_directory, "genomes"));
        writer_genome.save_genomes(info.genomes); 
      } 
    } 
  }


  INFO("MGRA log can be found here " << cfg::get().out_path_to_logger_file)
  INFO("Thank you for using MGRA!")

  return 0;
}
