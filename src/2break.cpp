#include "2break.h"

list< TwoBreak > TwoBreak::History;

// check if a 2-break creates a circular chromosome
bool TwoBreak::islinear(MBGraph& M) const {

    apply(M);

    for(int i = 0; i < 2; ++i) {
	vertex_t x;
	if( i==0 ) x = OldArc[0].first;
        else x = OldArc[0].second;

        if( x==Infty ) continue;

	for(auto ic = MultiColor.cbegin(); ic != MultiColor.cend(); ++ic) {

	    const partgraph_t& PG = M.get_local_graph(ic->first);
	    bool circular = false;

	    for(std::string y = M.get_adj_vertex(x); PG.defined(y);) {
		y = PG[y];

		if( y!=x ) y = M.get_adj_vertex(y);

		if( y==x ) {
		    circular = true;
                    break;
		}
	    }
		
	    if( !circular && PG.defined(x) ) {
		for(string y = x; PG.defined(y);) {
		    y = PG[y];

		    if(y!=x) y = M.get_adj_vertex(y);

		    if(y==x) {
			circular = true;
			break;
		    }
		}
	    }

	    if( circular ) {
		revert(M);
		return false;
	    }
	}
    }
    revert(M);
    return true;
}

bool TwoBreak::apply(MBGraph& M, bool record) const  {
   if (record) {
	outlog << " " << *this << " ";

#ifdef PREVENT_CIRCULAR
	if (!islinear(M)) {
           outlog << "REJECTED ";
	   return false;
	}
#endif

	History.push_back(*this);
   }
 
  for(auto ic = MultiColor.cbegin(); ic != MultiColor.cend(); ++ic) {
    for(size_t i = 0; i < 2; ++i) {
      if(OldArc[i].first != Infty && OldArc[i].second != Infty) {
	M.erase_edge(ic->first, OldArc[i].first, OldArc[i].second);
      } 
    }

      	    if(OldArc[0].first != Infty && OldArc[1].first != Infty) {
	        M.add_edge(ic->first, OldArc[0].first, OldArc[1].first);
	    }

      	    if(OldArc[0].second != Infty && OldArc[1].second != Infty) {
	        M.add_edge(ic->first, OldArc[0].second, OldArc[1].second);
	    }
	}
	return true;
}

void TwoBreak::applySingle(partgraph_t& SG) const {
	for(size_t i=0;i<2;++i) {
	    if(OldArc[i].first!=Infty && OldArc[i].second!=Infty) {
		SG.erase(OldArc[i].first,OldArc[i].second);
	    }
	}
	if(OldArc[0].first!=Infty && OldArc[1].first!=Infty) {
	    SG.insert(OldArc[0].first,OldArc[1].first);
	}
	if(OldArc[0].second!=Infty && OldArc[1].second!=Infty) {
	    SG.insert(OldArc[0].second,OldArc[1].second);
	}
    }

 void TwoBreak::normalize() {
	while( true ) {
	    if( OldArc[0].first > OldArc[0].second ) {
		OldArc[0] = make_pair( OldArc[0].second, OldArc[0].first );
		OldArc[1] = make_pair( OldArc[1].second, OldArc[1].first );
                continue;
	    }
	    if( OldArc[0].first > OldArc[1].first || OldArc[0].first > OldArc[1].second ) {
		arc_t temp = OldArc[0];
		OldArc[0] = OldArc[1];
		OldArc[1] = temp;
                continue;
	    }
	    break;
	}
    }

/*
If record == true we write elements in history
*/
void ApplyAll(MBGraph& M, const transform_t& T, bool record) {
    for(auto it = T.cbegin(); it != T.cend(); ++it) {
	it->apply(M, record);
    }
}

void RevertAll(MBGraph& M, const transform_t& T) {
    for(auto it = T.crbegin(); it != T.crend(); ++it) {
	it->revert(M);
    }
}

