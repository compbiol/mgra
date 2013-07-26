#include "mcolor.h"

bool Mcolor::includes(const Mcolor& second) const { 
	auto first1 = main_color.cbegin(); 
	auto first2 = second.cbegin(); 

	for (; first2 != second.cend(); ++first1) { 
		if (first1 == main_color.cend()) {
			return false;
		}

		if (first2->first < first1->first) { 
			return false;
		} 

		if (first1->first == first2->first) { 
			if  (first1->second < first2->second) {
				return false; 
			} else {
				++first2; 
			}
		} 
	} 
	return true;
} 

bool Mcolor::is_one_to_one_match() const { 
	for (auto it = main_color.cbegin(); it != main_color.cend(); ++it) { 
		if (it->second != 1) { 
			return false; 
		} 
	} 
		
	return true; 
} 

size_t Mcolor::how_much_includes(const Mcolor& second) const {
	size_t answer = 0;
	Mcolor current = *this; 
	while (current.includes(second)) {
		++answer; 
		current = Mcolor(current, second, Mcolor::Difference);
	}  	
	return answer; 
} 

void Mcolor::set_difference(const Mcolor& first, const Mcolor& second) { 
	auto first1 = first.cbegin(); 
	auto first2 = second.cbegin(); 
	auto result = std::inserter(main_color, main_color.begin());

	while (first1 != first.cend() && first2 != second.cend()) {
		if (first1->first == first2->first) { 
			if (first1->second > first2->second) { 
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
} 

