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

