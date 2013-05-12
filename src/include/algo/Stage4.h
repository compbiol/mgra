#ifndef STAGE_4_H_ 
#define STAGE_4_H_ 

inline bool stage4(MBGraph& graph) {

   outlog << "SplitBadColors is ON" << std::endl;
   graph.SplitBadColors = true; 

   Stage2 st2;
   bool isChanged = st2.stage2(graph);

   graph.SplitBadColors = false;
   outlog << "SplitBadColors is OFF" << std::endl;	 
	
   return isChanged;
} 

#endif
