#ifndef MCOLOR_H_
#define MCOLOR_H_

#include <iostream>
#include <algorithm> 
#include <map>
#include <set>
#include <utility>

#include "equivalence.h"

#define member(S,x) ((S).find(x)!=(S).end())

struct Mcolor {
	enum Construct {Difference, Union, Intersection};
	
	typedef std::map<size_t, size_t> my_map;
	typedef my_map::const_iterator citer;
	typedef my_map::iterator iter; 

	Mcolor() { 
	} 

	Mcolor(size_t i) {
		main_color.insert(std::make_pair(i, 1));
	}  

	Mcolor(const Mcolor& first, const Mcolor& second, const Construct& what) { 
		switch (what) { 
			case Difference: set_difference(first, second);	
					 break;
			case Union: set_union(first, second); 
		  		    	 break; 
			case Intersection: set_intersection(first, second);
					 break; 
		} 
	} 

	template<class colors_t>
	std::set<Mcolor> split_color(const colors_t& colors, bool split_bad_colors) const;
	
	bool is_good_multiedge() const;
	bool includes(const Mcolor& second) const;

	inline void insert(size_t i) {
		if (main_color.find(i) == main_color.end()) {  
			main_color.insert(std::make_pair(i, 1));
		} else {
			main_color[i] += 1;	
		} 
	} 

	inline size_t find(size_t i) { 
		return main_color[i];					
	} 

	inline bool mymember(size_t j) const { 
		return (main_color.find(j) != main_color.end());
	} 

	inline bool empty() const { 
		return main_color.empty();
	} 

	inline size_t size() const { 
		return main_color.size();
	} 
	
	inline iter begin() { 
		return main_color.begin();
	} 

	inline iter end() { 
		return main_color.end();
	} 

	inline citer cbegin() const { 
		return main_color.cbegin();
	} 

	inline citer cend() const { 
		return main_color.cend();
	} 

	bool operator > (const Mcolor& C) const { 
		return (main_color > C.main_color);
	} 

	bool operator < (const Mcolor& C) const {
		return (main_color < C.main_color); 
	} 

	bool operator == (const Mcolor& C) const { 
		return (main_color == C.main_color);	
	} 

	bool operator != (const Mcolor& C) const { 
		return (main_color != C.main_color); 
	}
	
private:
	void set_difference(const Mcolor& first, const Mcolor& second);
	void set_union(const Mcolor& first, const Mcolor& second);
	void set_intersection(const Mcolor& first, const Mcolor& second); 
private: 
	my_map main_color;
};

/*
SplitColor(Q) представляет Q в виде дизъюнктного объединения T-consistent мультицветов, т.е. Q = Q1 U ... U Qm
где каждый Qi является T-consistent и все они попарно не пересекаются. SplitColor(Q) возвращает множество { Q1, Q2, ..., Qm }
(в частности, когда Q является T-consistent, имеем m=1 и Q1=Q).
Теперь, когда SplitBadColors = true, то и ребро (x,y) имеет мультицвет Q, то MBG.mulcols(x) будет содежать вместо (Q,x) пары:
(Q1,y), (Q2,y), ..., (Qm,y)
*/
template<class colors_t>
std::set<Mcolor> Mcolor::split_color(const colors_t& colors, bool split_bad_colors) const {
    std::set<Mcolor> S;

    if (colors.is_vec_T_color(*this)) {
	S.insert(*this);
	return S;
    }

    if (!split_bad_colors) {
        return S;
    }

    equivalence<size_t> EQ;
    for(auto iq = this->cbegin(); iq != this->cend(); ++iq) { 
	EQ.addrel(iq->first, iq->first);
    } 

    for(auto ic = colors.cbegin_T_color(); ic != colors.cend_T_color(); ++ic) {
	Mcolor C(*ic, *this, Mcolor::Intersection);
	if (C.size() >= 2 && C.size() == ic->size() ) {
	    for (auto iq = C.begin(); iq != C.end(); ++iq) {
		EQ.addrel(iq->first, C.begin()->first);
	    }
	}
    }

    EQ.update();
    std::map<size_t, Mcolor> cls; 
    EQ.get_eclasses(cls);
    for(auto ic = cls.cbegin(); ic != cls.cend(); ++ic) {
	S.insert(ic->second);
    }
    return S;
} 

#endif
