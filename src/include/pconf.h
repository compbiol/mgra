/* 
** Module: MGRA Configurations support
** Version: 1.1
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008,12 by Max Alekseyev <maxal@cse.sc.edu> 
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/

#ifndef PCONF_H_
#define PCONF_H_

#include <unordered_map>
#include <vector>
#include <list>
#include <sstream>
#include <iostream>
#include <cassert>

#include <mcolor.h>

template <typename T>
__attribute__((always_inline)) inline std::string toString(T val) { 
	std::ostringstream oss; 
	oss << val; 
	return oss.str();
}  

__attribute__((always_inline)) inline bool check_symbol(char c) { 
   if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9')) { 
     return true;
   } else { 
     return false;
   } 
} 


/*This structures containes information from *.cfg file */
struct ProblemInstance {
	ProblemInstance(const std::unordered_map<std::string, std::vector<std::string> >& input); 

	inline std::vector<std::string>::const_iterator cbegin_name_genome() const {  ///FIXME: USING in Wstats.txt
		return priority_name.cbegin();
	} 

	inline std::vector<std::string>::const_iterator cend_name_genome() const { 
		return priority_name.cend();
	} 

	inline size_t get_count_genomes() const { 
		return priority_name.size();
	} 

	inline static std::string get_name(size_t i) { 
		return priority_name[i];
	} 

	inline static bool member_name(const std::string& i) { 
		return (genome_number.find(i) != genome_number.end());
	} 

	inline static size_t get_number(const std::string& str) {
		assert (member_name(str)); 
		return genome_number.find(str)->second;
	} 
	
	inline std::string get_blk_format() const { 
		return block_format;
	} 

	inline std::string get_blk_file() const { 
		return block_file;		
	} 

	inline std::string get_graphname() const { 
		return graphfname;		
	} 

	inline std::string get_colorscheme() const { 
		return colorscheme;		
	} 

	inline std::vector<std::string>/*std::vector<BinaryTree<std::string> >*/ get_trees() const { 
		return trees;
	} 

	inline size_t get_stages() const { 
		return stages;
	} 

	inline Mcolor get_target() const { 
		return target;		
	} 

	//inline Mcolor get_newtarget() const { 
	//	return color_target;		
	//} 

	inline std::list<std::vector<std::string> > get_completion() const { //del and write convert two break
		return completion;
	} 

	inline size_t RGB_color_size() { 
		return RGBcolors.size();
	} 

	inline std::string get_RGBcolor(int index) const { 
		return RGBcolors[index];
	} 

	inline int get_RGBcoeff() const { 
		return RGBcoeff;
	} 
private: 
	void init_basic_rgb_colors(bool flag = true); 
private:
	static std::vector<std::string> priority_name;
	//std::vector<std::string> priority_name;  
	static std::unordered_map<std::string, size_t> genome_number;
	//std::unordered_map<std::string, size_t> genome_number;
	static std::unordered_map<size_t, std::string> number_genome;
	//std::unordered_map<size_t, std::string> number_genome;  

	std::string block_format; 	
	std::string block_file;

	std::string graphfname;
	std::string colorscheme;

	//std::vector<BinaryTree<std::string> > trees;
	std::vector<std::string> trees;

	size_t stages;

	Mcolor target; 		
	//std::string target;

	std::list<std::vector<std::string> > completion;

	/*about color refactor in future. What*/
	std::vector<std::string> RGBcolors;
	int RGBcoeff; 
};

#endif
