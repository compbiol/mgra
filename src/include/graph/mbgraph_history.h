#ifndef MBGRAPH_HISTORY_H_
#define MBGRAPH_HISTORY_H_

#include "mbgraph_colors.h"
#include "2break.h"

#define member(S,x) ((S).find(x)!=(S).end())

template<class mcolor_t>
struct mbgraph_with_history: public mbgraph_with_colors<mcolor_t> { 

	mbgraph_with_history(const std::vector<Genome>& genomes, const ProblemInstance& cfg)
	: mbgraph_with_colors<mcolor_t>(genomes, cfg) 
	{ 
	} 

	inline void insert_twobreak(const TwoBreak<mbgraph_with_history, mcolor_t>& break2) {
	 	history.push_back(break2);
	} 

	inline typename std::list<TwoBreak<mbgraph_with_history, mcolor_t> > get_history() const { //FIXME: DEL
		return history; 
	} 

	inline typename std::list<TwoBreak<mbgraph_with_history, mcolor_t> >::const_iterator begin_history() const { 
		return history.cbegin();
	} 
	
	inline typename std::list<TwoBreak<mbgraph_with_history, mcolor_t> >::const_iterator end_history() const { 
		return history.cend(); 
	} 
private: 
	std::list<TwoBreak<mbgraph_with_history, mcolor_t> > history;	
};

typedef std::list<TwoBreak<mbgraph_with_history<Mcolor>, Mcolor> > transform_t;

#endif
