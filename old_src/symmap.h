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

#include <map>
#include <functional>
#include <iostream>

using namespace std;

template<class item_class, class order = less<item_class> >
    class symmap : public map<item_class,item_class,order>
{
    int card;

    typedef map<item_class,item_class,order> B;

public:

    typedef typename B::const_iterator const_iterator;

    using B::end;
    using B::find;

    symmap() {
        card = 0;
    }

    bool defined(const item_class& x) const {
        return find(x) != end();
    }

    void insert(const item_class& x,const item_class& y) {
        if( defined(x) || defined(y) ) {
            cerr << endl << "symmap::insert() error: redefining ";
            if( defined(x) ) cerr << "(" << x << "," << B::operator[](x) << ") ";
	    if( defined(y) ) cerr << "(" << y << "," << B::operator[](y) << ") ";
            cerr << "with (" << x << "," << y << ")" << endl;
	    abort();
        }
        B::operator[](x) = y;
        B::operator[](y) = x;
        card++;
    }

    void erase(const item_class& x) {
        if(!defined(x)) {
            cerr << endl << "symmap::erase() error: unmapped element (" << x << ")" << endl;
            abort();
        }
        B::erase(B::operator[](x));
        B::erase(x);
        card--;
    }

    const item_class& operator [] (const item_class& x) const {
        const const_iterator ix = find(x);
	if(ix == end()) {
            cerr << "symmap::operator[] error: undefined element " << x << endl;
            abort();
        }
        return ix->second;
    }

    void erase(const item_class& x,const item_class& y) {
        if(!defined(x) || B::operator[](x)!=y) {
            cerr << "symmap::erase() error: unmapped pair (" << x << "," << y << ")" << endl;
            abort();
        }
        B::erase(x);
        B::erase(y);
        card--;
    }

    int size() const {
        return card;
    }

    bool empty() const {
        return card==0;
    }

    void clear() {
        B::clear();
        card = 0;
    }

    void pop(item_class& x,item_class& y) {
        if(!card) {
            cerr << "symmap::extract error: empty symmap!\n";
            abort();
        }
        x = B::begin()->first;
        y = B::begin()->second;
        B::erase(x);
        B::erase(y);
        card--;
    }

/*
    bool istrue(const item_class& x,const item_class& y) const {
        if(!defined(x)) return false;
        const_iterator ci = find(x);
        neighborhood::const_iterator t = ci->second.find(y);
        return t != ci->second.end();
    }

    pair<item_class,item_class> pop() {
        item_class x = begin()->first, y = *(begin()->second.begin()); //y = *(operator[](x).begin());
        erase(x,y);
        return make_pair(x,y);
    }

    set<item_class> adjacent(const item_class& x) const {
        set<item_class> res;
        if(defined(x)) {
            const_iterator ci = find(x);
            res = ci->second;
        }
        return res;
    }
*/
};

#endif
