#include "mcolor.h"

bool Mcolor::includes(const Mcolor& second) const { 
	auto first1 = main_color.cbegin(); 
	auto first2 = second.cbegin(); 

	for (; first2 != second.cend(); ++first1) { 
		if (first1 == main_color.cend()) {
			return false;
		}

		if ((first1->first == first2->first) && (first1->second < first2->second)) { 
			return false; 
		} 
 
		if (*first2 < *first1) { 
			return false;
		} 

		if (*first1 == *first2) {
			++first2;
		} 
	} 
	return true;
	//return std::includes(main_color.cbegin(), main_color.cend(), second.main_color.cbegin(), second.main_color.cend());
} 

bool Mcolor::is_good_multiedge() const { 
	for (auto it = main_color.cbegin(); it != main_color.cend(); ++it) { 
		if (it->second != 1) { 
			return false; 
		} 
	} 
		
	return true; 
} 

/*
SplitColor(Q) представляет Q в виде дизъюнктного объединения T-consistent мультицветов, т.е. Q = Q1 U ... U Qm
где каждый Qi является T-consistent и все они попарно не пересекаются. SplitColor(Q) возвращает множество { Q1, Q2, ..., Qm }
(в частности, когда Q является T-consistent, имеем m=1 и Q1=Q).
Теперь, когда SplitBadColors = true, то и ребро (x,y) имеет мультицвет Q, то MBG.mulcols(x) будет содежать вместо (Q,x) пары:
(Q1,y), (Q2,y), ..., (Qm,y)
*/
/*std::set<Mcolor> Mcolor::split_color(const ColorsGraph<Mcolor>& colors, bool split_bad_colors = false) const {
    std::set<Mcolor> S;

    if (member(colors.DiColor, *this)) {
	S.insert(*this);
	return S;
    }

    if (!split_bad_colors) {
        return S;
    }

    equivalence <size_t> EQ;
    for(auto iq = this->cbegin(); iq != this->cend(); ++iq) { 
	EQ.addrel(iq->first, iq->first);
    } 

    for(auto ic = colors.DiColor.begin(); ic != colors.DiColor.end(); ++ic) {
	Mcolor C(*ic, *this, Mcolor::Intersection);
	if (C.size() >= 2 && C.size() == ic->size() ) {
	    for (auto iq = C.begin(); iq != C.end(); ++iq) {
		EQ.addrel(iq->first, C.begin()->first);
	    }
	}
    }

    EQ.update();
    std::map <size_t, Mcolor > cls; //FIXME because cls - is Mcolor
    EQ.get_eclasses(cls);
    for(auto ic = cls.cbegin(); ic != cls.cend(); ++ic) {
	S.insert(ic->second);
    }
    return S;
}*/

void Mcolor::set_difference(const Mcolor& first, const Mcolor& second) { 
	auto first1 = first.cbegin(); 
	auto first2 = second.cbegin(); 
	auto result = std::inserter(main_color, main_color.begin());

	while (first1 != first.cend() && first2 != second.cend()) {
		if (first1->first == first2->first) { 
			if (first1->second - first2->second != 0) { 
				*result = std::make_pair(first1->first, first1->second - first2->second);
			} 
			++first1;
			++first2;
		} else if (*first1 < *first2) {
			*result = *first1; 
			++result; 
			++first1; 
		} else if (*first2 < *first1) { 
			++first2;
		} 
	}

	for (; first1 != first.cend(); ++first1, ++result) { 
		*result = *first1;
	} 
	//std::set_difference(first.cbegin(), first.cend(), second.cbegin(), second.cend(), std::inserter(main_color, main_color.begin()));
} 

void Mcolor::set_union(const Mcolor& first, const Mcolor& second) { 
	auto first1 = first.cbegin(); 
	auto first2 = second.cbegin(); 
	auto result = std::inserter(main_color, main_color.begin());

	while (true) { 
		if (first1 == first.cend()) { 
			for (; first2 != second.cend(); ++first2, ++result) { 
				*result = *first2;
			} 
			break; 
		} 

		if (first2 == second.cend()) { 	
			for (; first1 != first.cend(); ++first1, ++result) { 
				*result = *first1;
			} 		
			break;
		} 

		if (first1->first == first2->first) { 
			*result = std::make_pair(first1->first, first1->second + first2->second);	
			++first1; 
			++first2; 
		} else if (*first1 < *first2) {
			*result = *first1; 
			++first1; 
		} else if (*first2 < *first1) { 
			*result = *first2; 
			++first2;			
		}  
		++result;
	} 	
	//std::set_union(first.cbegin(), first.cend(), second.cbegin(), second.cend(), std::inserter(main_color, main_color.begin()));
} 

void Mcolor::set_intersection(const Mcolor& first, const Mcolor& second) { 
	auto first1 = first.cbegin(); 
	auto first2 = second.cbegin(); 
	auto result = std::inserter(main_color, main_color.begin());

	while (first1 != first.cend() && first2 != second.cend()) {
		if (first1->first == first2->first) { 
			*result = std::make_pair(first1->first, std::min(first1->second, first2->second)); 		
			++result; 
			++first1;
			++first2;
		} else if (*first1 < *first2) { 
			++first1;
		} else if (*first2 < *first1) { 
			++first2;
		} 
	}
	//std::set_intersection(first.cbegin(), first.cend(), second.cbegin(), second.cend(), std::inserter(main_color, main_color.begin()));
} 

