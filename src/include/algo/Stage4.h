#ifndef STAGE_4_H_ 
#define STAGE_4_H_ 

inline bool stage4(MBGraph& graph, bool canformQ) {

   outlog << "SplitBadColors is ON" << std::endl;
   graph.SplitBadColors = true; 

   Stage2<MBGraph> st2(graph, canformQ);
   bool isChanged = st2.stage2();

   graph = st2.get_graph();
   graph.SplitBadColors = false;
   outlog << "SplitBadColors is OFF" << std::endl;	 
	
   return isChanged;
} 

#endif
