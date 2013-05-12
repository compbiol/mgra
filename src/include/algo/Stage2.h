#ifndef STAGE_2_H_ 
#define STAGE_2_H_ 

#include <string>

#include "mpbgraph.h"
#include "2break.h"

extern bool canformQoo; // safe choice, at later stages may change to false
bool canformQ(const std::string& x, const Mcolor& Q);

struct Stage2 { 
	/*Stage2(bool canform)
	:canformQoo(canform) { 
	} */

	bool stage2(MBGraph& graph);

	/*
	canformQ(x,Q) говорит, можно ли из мультицветов мультиребер инцидентных вершине x образовать мультицвет Q.
	can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
	if return false, then Q cannot be formed
	if true - who knows...
	*/
	//bool canformQ(const std::string& x, const Mcolor& Q);


	//bool canformQoo  = true; 		 
}; 	

#endif
