#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "genome_match.h" //FIXME REMOVE LATER

#include "graph/breakpoint_graph.hpp"
#include "writer/Wstats.h"
#include "writer/Wdots.h"

template<class graph_t>
struct Algorithm { 
  Algorithm(std::shared_ptr<graph_t> const & gr, ProblemInstance<typename graph_t::mcolor_type> const & cfg) 
  : graph(gr) 
  , rounds(cfg.get_max_number_of_split())
  , max_size_component(cfg.get_size_component_in_brutforce())
  , stages(cfg.get_stages())
  , m_completion(cfg.get_completion())
  , write_dots(cfg)
  {
  } 

  void init_writers(fs::path const & work_dir, std::string const & colorscheme, std::string const & graphname, bool debug) { 
    write_stats.open(work_dir, "stats.txt");
    write_dots.init(work_dir, colorscheme, graphname, debug);
  }

  void convert_to_identity_bgraph();
  
private: 
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;
  
  struct Stage {
    explicit Stage(std::shared_ptr<graph_t> const & gr) 
    : graph(gr)
    {
    }

    virtual bool do_action() = 0;
    virtual std::string get_name() = 0;
    virtual ~Stage() 
    {
    }
  protected:
    std::shared_ptr<graph_t> graph;
  };

public:
  struct Balance; 	

private:
  struct ProcessSimplePath;
  
  struct ProcessFairEdges; 

  struct ProcessClone;

  struct IncreaseNumberComponents;

  struct BruteForce;
  
  //Stage 6: brutoforce stage
  bool stage6();
  
  std::multimap<size_t, arc_t> create_minimal_matching(std::set<vertex_t> const & vertex_set); 
  std::multimap<size_t, arc_t> wrapper_create_minimal_matching(std::set<vertex_t> const & vertex_set);

  size_t take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y);
  
  /*EXPERIMENTAL STAGE OR EXISTS BUT NOT USED */
  //Stage -: process tandem duplication events
  bool stage5_4();
  //bool stage71();
  //stage5_3()    
  //stage4_td(); 
  
private: 
  std::shared_ptr<graph_t> graph; 

  size_t const rounds;
  size_t const max_size_component; 
  size_t const stages;
  std::list<twobreak_t> const m_completion;

  writer::Wstats write_stats;
  writer::Wdots<graph_t, ProblemInstance<mcolor_t> > write_dots;
};

template<class graph_t>
void Algorithm<graph_t>::convert_to_identity_bgraph() {
  std::unordered_set<size_t> print_dots;
  size_t stage = 0; 
  bool isChanged = false;
  bool process_compl = true; 

  std::vector<std::shared_ptr<Stage> > algorithm(5); 
  algorithm[0] = std::shared_ptr<Stage>(new Balance(graph));
  algorithm[1] = std::shared_ptr<Stage>(new ProcessSimplePath(graph));
  algorithm[2] = std::shared_ptr<Stage>(new ProcessFairEdges(graph));
  algorithm[3] = std::shared_ptr<Stage>(new ProcessClone(graph));
  algorithm[4] = std::shared_ptr<Stage>(new IncreaseNumberComponents(graph));
  
  auto const saveInfoLambda = [&](size_t st) -> void { 
    if ((print_dots.count(st) == 0) && !isChanged) {
      print_dots.insert(st);
      Statistics<graph_t> stat(graph);
      write_stats.print_all_statistics(st, stat, *graph);
      write_dots.save_dot(*graph, st);
    } 
  };

  saveInfoLambda(stage++);
  write_dots.write_legend_dot();

  if (stages >= 1) { 
    graph->update_number_of_splits(rounds);  
    std::cerr << "Stage " << stage << ": " << algorithm[0]->get_name() << std::endl; 
    isChanged = algorithm[0]->do_action();
    saveInfoLambda(stage++);
  }

  isChanged = true;
  while(isChanged) {
    isChanged = false; 
    stage = 2;

    for (size_t i = 1; i <= rounds && !isChanged; ++i) {   
      std::cerr << "Rounds " << i << std::endl;

      graph->update_number_of_splits(i);

      for (size_t j = 1; j < 5 && !isChanged; ++j) {
        std::cerr << "Stage " << stage << ": " << algorithm[j]->get_name() << std::endl;
        write_dots.save_dot(*graph, 100);      
        isChanged = algorithm[j]->do_action();
        saveInfoLambda(stage++);
      } 
    } 

    if (graph->get_canformQoo() && !isChanged) { 
      std::cerr << "Change canformQoo "  << std::endl;
      graph->set_canformQoo(false);
      isChanged = true;
    }

    /*if (is_mobile_irregular_edge && !isChanged) { 
      std::cerr << "Change is process irregular edge " << std::endl;
      canformQoo = true; // more flexible
      is_mobile_irregular_edge = false; //more flexible 
      isChanged = true;
    } */
    
    /*if (!isChanged) {
      std::cerr << "Stage 5: Experement stage" << stage << std::endl;
      graph->update_number_of_splits(3);
      isChanged = stage5_4();
      saveInfoLambda(stage++);
    }*/

    if (process_compl && !m_completion.empty() && !isChanged) {     
      //std::cerr << "Manual Completion Stage" << std::endl;
      for(auto il = m_completion.cbegin(); il != m_completion.cend(); ++il) {
        graph->apply(*il);
      }

      process_compl = false;
      isChanged = true;
    }  

    if ((max_size_component != 0) && !isChanged) {
      std::cerr << "Brute force stage: " << std::endl;
      graph->update_number_of_splits(3);
      isChanged = stage6();
      saveInfoLambda(stage++);
    }
  }	
  
  //graph->print();
  graph->update_number_of_splits(3);
  write_dots.save_final_dot(*graph);
}          

#include "algo/Stage1.h" 
#include "algo/Stage2.h"
#include "algo/Stage3.h"
//#include "algo/Stage4.h"
#include "algo/Stage5.h"
#include "algo/Stage6.h"
#include "algo/Stage7.h"
#include "algo/ExperemStages.h"


#endif
