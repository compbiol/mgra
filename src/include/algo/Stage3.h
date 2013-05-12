#ifndef STAGE_3_H_
#define STAGE_3_H_

//#include <unredered_set>
#include <string>

#include "mpbgraph.h"
#include "2break.h"

bool stage3(MBGraph& graph); 

// cut the graph into connected components
bool stage3_1(MBGraph& graph); 

 // search for 4-cycles
bool stage3_2(MBGraph& graph);

#endif
