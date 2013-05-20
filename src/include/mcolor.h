#ifndef MCOLOR_H_
#define MCOLOR_H_

#include <iostream>
#include <algorithm> 
#include <map>
#include <set>
#include <utility>

//#include "graph_colors.h"
//#include "utility/equivalence.h"

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

	//std::set<Mcolor> split_color(ColorsGraph<Mcolor>& colors, bool split_bad_colors) const;	

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



#endif
