/*
** Module: Symmetric maps ver. 1.5
**
** This file is part of the
** Multiple Genome Rearrangements and Ancestors (MGRA)
** reconstruction software.
**
** Copyright (C) 2008,09 by Max Alekseyev <maxal@cse.sc.edu>
** 
** This program is free software; you can redistribute it and/or
** modify it under the terms of the GNU General Public License
** as published by the Free Software Foundation; either version 2
** of the License, or (at your option) any later version.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html
*/

#ifndef SYMMAP_H
#define SYMMAP_H

#include <unordered_map>
#include <functional>
#include <iostream>

template<class item_class, class Hash = hash<item_class> >
struct symmap {
	private:
		int card;
		std::unordered_map<item_class, item_class, Hash > main_map;
	public:
		symmap() {
			card = 0;
		}

		inline bool defined(const item_class& x) const {
			return main_map.find(x) != end();
		}

		void insert(const item_class& x, const item_class& y) {
			if (defined(x) || defined(y)) {
				std::cerr << std::endl << "symmap::insert() error: redefining ";
				if (defined(x)) { 
					std::cerr << "(" << x << "," << main_map.find(x) << ") ";
				} 
				if (defined(y)) { 
					std::cerr << "(" << y << "," << main_map.find(y) << ") ";
				} 
				std::cerr << "with (" << x << "," << y << ")" << std::endl;
				abort();
		        }
			main_map.insert(std::make_pair(x, y));
			main_map.insert(std::make_pair(y, x));
        		++card;
		}

		void erase(const item_class& x) {
			if (!defined(x)) {
				std::cerr << std::endl << "symmap::erase() error: unmapped element (" << x << ")" << std::endl;
				abort();
			}
			main_map.erase(main_map.find(x));
			main_map.erase(x);
			--card;
		}

		const item_class& operator[] (const item_class& x) const {
			auto ix = main_map.find(x);
			if (ix == main_map.end()) {
				std::cerr << "symmap::operator[] error: undefined element " << x << std::endl;
				abort();
			}
			return ix->second;
		}

		void erase(const item_class& x,const item_class& y) {
			if (!defined(x) || main_map.find(x) != y) {
				std::cerr << "symmap::erase() error: unmapped pair (" << x << "," << y << ")" << std::endl;
				abort();
			}
			main_map.erase(x);
			main_map.erase(y);
			--card;
		}

		inline int size() const {
			return card;
		}

		inline bool empty() const {
			return (card == 0);
		}

		void clear() {
			main_map.clear();
			card = 0;
		}

		void pop (item_class& x, item_class& y) {
			if (card == 0) {
				std::cerr << "symmap::extract error: empty symmap!" << std::endl;
				abort();
			}
			x = main_map.begin()->first;
			y = main_map.begin()->second;
			main_map.erase(x);
			main_map.erase(y);
			--card;
		}
};

#endif
