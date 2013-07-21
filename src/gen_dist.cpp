#include "gen_dist.h"

std::map<edge_t, double> BEC; // breakpoint reuse on black edges
std::map<edge_t, double> GEC; // breakpoint reuse on gray edges, 
std::map<orf_t, double> BPR; // and endpoints
std::set<orf_t> USE; //used verticed 
std::set<orf_t> USE2; //and used at least 2 times


// returns: # trivial cycles, # (true) synteny blocks, d2, d3
std::vector<size_t> genome_dist(const partgraph_t& BE, const partgraph_t& GE, const partgraph_t& OE, const bool closed) {
	std::vector<size_t> res(4);

BEC.clear(); 
GEC.clear();
BPR.clear();
USE.clear();
USE2.clear();

/*
    L.clear();
    Lbb.clear();
    Lgg.clear();
    Lbg.clear();
*/
    map<int,int> L, Lbb, Lgg, Lbg;
    int  ncyc, bb, gg, bg;
    {   // count chains
        ncyc = bb = gg = bg = 0;
	set<string> processed;
	for(partgraph_t::const_iterator oi=OE.begin();oi!=OE.end();++oi) for(int j=0;j<2;++j) {
	    const string& gene = (j==0)?(oi->first):(oi->second);

	    if(processed.find(gene)!=processed.end()) continue;
	    if(!BE.defined(gene) && !GE.defined(gene)) {
		processed.insert(gene); 
/*		
		BPR[gene] += 2; USE.insert( gene );
*/
		//cout << "Isolated " << gene << endl;
		Lbg[1]++;
		bg++;
		continue;
	    }

//set<orf_t> CurrG;
	    int len = 0;
	    if(!BE.defined(gene) && GE.defined(gene)) {
		string y = gene;
		while(processed.find(y)==processed.end()) {
		    processed.insert(y); //CurrG.insert(y);
		    len++;
		    if (GE.defined(y)) y = GE[y];
		    else {
			bg++;
			// cout << "bg:" << len << " ";
			Lbg[len]++;
			break;
		    }
		    processed.insert(y); //CurrG.insert(y);
		    len++;
		    if(BE.defined(y)) y = BE[y];
		    else {
			gg++;
			// cout << "gg:" << len << " ";
			Lgg[len]++;
			break;
		    }
		}
	    }
	    if(BE.defined(gene) && !GE.defined(gene)) {
		string y = gene;
		while(processed.find(y)==processed.end()) {
		    processed.insert(y); //CurrG.insert(y);
		    len++;
		    if(BE.defined(y)) y = BE[y];
		    else {
			bg++;
			// cout << "bg:" << len << " ";
			Lbg[len]++;
			break;
		    }
		    processed.insert(y); //CurrG.insert(y);
		    len++;
		    if(GE.defined(y)) y = GE[y];
		    else {
			bb++;
			// cout << "bb:" << len << " ";
			Lbb[len]++;
			break;
		    }
		}
	    }
/*
	    for(set<orf_t>::const_iterator ic=CurrG.begin();ic!=CurrG.end();++ic) {
		BPR[*ic] += (len-4)/(2.*len) + 2;
		USE.insert( *ic );
	    }	    
*/
	}		    
#if 1

	for(partgraph_t::const_iterator oi=OE.begin();oi!=OE.end();++oi) for(int j=0;j<2;++j) {
	    string y = (j==0)?(oi->first):(oi->second);

	    if(processed.find(y)!=processed.end()) continue;
	    ncyc++;
	    int len = 0;

    set< edge_t > BES, GES;

    	    while(processed.find(y)==processed.end()) {
		
    edge_t e1; e1.insert(y); e1.insert(BE[y]);
    BES.insert(e1);
		
	        processed.insert(y);
		y = BE[y];
	
    edge_t e2; e2.insert(y); e2.insert(GE[y]);
    GES.insert(e2);

		processed.insert(y);
		y = GE[y];
		len+=2;
	    }
	    //cout << "cyc:" << len << " ";
	    L[len]++;

	      for(set<edge_t>::const_iterator is=BES.begin();is!=BES.end();++is) {
		if( len>4 ) {
	    	  BEC[*is] = (len-4)/(double)len;
		  BPR[*(is->begin())] += (len-4)/(2.*len);
		  BPR[*(is->rbegin())] += (len-4)/(2.*len);
		}
		if( len>=4 ) {
		    USE.insert( *(is->begin()) );
		    USE.insert( *(is->rbegin()) );
		}    
		if( len>4 ) {
		    USE2.insert( *(is->begin()) );
		    USE2.insert( *(is->rbegin()) );
		}    
	      }
	      for(set<edge_t>::const_iterator is=GES.begin();is!=GES.end();++is) {
		if( len>4 ) {
	          GEC[*is] = (len-4)/(double)len;
		  BPR[*(is->begin())] += (len-4)/(2.*len);
		  BPR[*(is->rbegin())] += (len-4)/(2.*len);
	        }		
		if( len>=4 ) {
		    USE.insert( *(is->begin()) );
		    USE.insert( *(is->rbegin()) );
		}    
		if( len>4 ) {
		    USE2.insert( *(is->begin()) );
		    USE2.insert( *(is->rbegin()) );
		}    
	      }


	}
	//cout << endl;
//	cout << "# bg/bb/gg: " << bg << " " << bb << " " << gg << endl;
//	cout << "# cycles: " << ncyc << "\ttrivial: " << L[2] << endl;

res[0] = L[2];
    	
	int nblk = OE.size() - L[2];
	ncyc -= L[2];
	L.erase(2);

res[1] = nblk;
	
//	cout << "# blocks: " << nblk << endl;
//	cout << endl;
	
	int tot = 0, tL = 0, ncyc_odd = 0;

//	cout << "cycles: ";
	for(map<int,int>::const_iterator il=L.begin();il!=L.end();++il) {
//	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    tL += il->second;
	    if((il->first/2)%2) ncyc_odd += il->second;
	}
//	cout << "(total vertices: " << tot << " " << "cycles: " << tL << ")" << endl;
//	cout << "bodd/godd: " << ncyc_odd << "\tbeven/geven: " << tL - ncyc_odd << endl;
//	cout << endl;
if( closed ) {
    res[2] = nblk - tL;
    res[3] = (nblk - ncyc_odd) / 2;
}

if( !closed ) {
//	cout << "bb-chains: ";
	tot = 0; 
	int nlbb = 0, nlbb_bodd = 0;
	for(map<int,int>::const_iterator il=Lbb.begin();il!=Lbb.end();++il) {
//	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlbb += il->second;
	    if( (il->first/2)%2 ) nlbb_bodd += il->second;
	}
//	cout << "(total vertices: " << tot << " " << "chains: " << nlbb << ")" << endl;
//	cout << "bodd/geven: " << nlbb_bodd << "\tbeven/godd: " << nlbb - nlbb_bodd << endl;
//	cout << endl;
	
//	cout << "gg-chains: ";
	tot = 0;
	int nlgg = 0, nlgg_bodd = 0;
	for(map<int,int>::const_iterator il=Lgg.begin();il!=Lgg.end();++il) {
//	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlgg += il->second;
	    if( !((il->first/2)%2) ) nlgg_bodd += il->second;
	}
//	cout << "(total vertices: " << tot << " " << "chains: " << nlgg << ")" << endl;
//	cout << "bodd/geven: " << nlgg_bodd << "\tbeven/godd: " << nlgg - nlgg_bodd << endl;
//	cout << endl;

//	cout << "bg-chains: ";
	tot = 0;
	int nlbg = 0, nlbg_bodd = 0;
	for(map<int,int>::const_iterator il=Lbg.begin();il!=Lbg.end();++il) {
//	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlbg += il->second;
	    if( ((il->first-1)/2)%2 ) nlbg_bodd += il->second;
	}
//	cout << "(total vertices: " << tot << " " << "chains: " << nlbg << ")" << endl;
//	cout << "bodd/godd: " << nlbg_bodd << "\tbeven/geven: " << nlbg - nlbg_bodd << endl;
//	cout << endl;
	
	int d1 = (nlgg - nlgg_bodd) - (2*nlbg_bodd - nlbg)/2;
	int d2 = nlbb_bodd - (2*nlbg_bodd - nlbg)/2;
	if( d1>0 ) d1 %= 2; else d1 = 0;
	if( d2>0 ) d2 %= 2; else d2 = 0;

//	cout << "B_2: " << nblk - ncyc - max(1,nlbg/2) - nlbb << "\t" << nblk - ncyc - max(1,nlbg/2) - nlgg << endl;

	res[2] = max( nblk - ncyc - max(1,nlbg/2) - nlbb, nblk - ncyc - max(1,nlbg/2) - nlgg );

	int L3PQ = ncyc_odd + nlbb_bodd + d1 + max(0,abs(2*nlbg_bodd - nlbg)/2 - (nlgg - nlgg_bodd));
	int L3QP = ncyc_odd + (nlgg - nlgg_bodd) + d2 + max(0,abs(2*nlbg_bodd - nlbg)/2 - nlbb_bodd);
	
//	cout << "B_3: " << (nblk - L3PQ)/2 << "\t" << (nblk - L3QP)/2 << endl;

	res[3] = max( (nblk - L3PQ)/2, (nblk - L3QP)/2 );

//	cout << "B_23: " << (3*nblk - 2*nlbb - 2*ncyc - nlbg - L3PQ)/2 << "\t" 
//	<< (3*nblk - 2*nlgg - 2*ncyc - nlbg - L3QP)/2 << endl;
//	cout << endl;
}	
//	cout << endl;
#endif
    }		    
    return res;
}
