/* 
** Module: Equivalence relation ver. 1.3
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008,09 by Max Alekseyev <maxal@cse.sc.edu> 
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/

#ifndef EQUIV_H
#define EQUIV_H

#include<map>
#include<functional>

template< class Item, class Cmp=std::less<Item> >
    class equivalence : public map<Item,Item,Cmp> {

    typedef map<Item,Item,Cmp> map_t;

    static Cmp cmp;

public:

    /*
    typedef typename map_t::iterator iterator;
    typedef typename map_t::const_iterator const_iterator;
    typedef typename map_t::reverse_iterator reverse_iterator;
    typedef typename map_t::const_reverse_iterator const_reverse_iterator;
    */

    //typedef map<Item, set<Item>, Cmp> eclass_t;

    using map_t::begin;
    using map_t::end;
    using map_t::find;

    //typedef set<Item> class_t;


    // creates equivalence relation
    // equivalence();
    // destructor
    //~equivalence();


    // introduce a relation between two specified integers
    void addrel(const Item&,const Item&);
    void addrel(const pair<Item,Item>& p) { addrel(p.first,p.second); }

    // return `true' iff two integers are equivalent
    bool isequiv(const Item&,const Item&);
    bool isequiv(const pair<Item,Item>& p) { return isequiv(p.first,p.second); }

    const Item& operator [] (const Item&);

    void insert(const Item& x) {
	operator[](x);
    }

    bool defined(const Item&) const;

    void update ();

    int classes();

    //set<int> equivclass(int);

    /* !! get_eclasses() const must work on updated class ONLY !! */
    //template<class eclass_t>
    //void get_eclasses(map<Item,eclass_t,Cmp>&) const;

    template<class eclass_t>
    void get_eclasses(map<Item,eclass_t,Cmp>&);

};

template<class Item, class Cmp>
Cmp equivalence<Item,Cmp>::cmp;

template<class Item, class Cmp>
bool equivalence<Item,Cmp>::defined(const Item& x) const {
    return find(x) != end();
}

template<class Item, class Cmp>
const Item& equivalence<Item,Cmp>::operator [] (const Item& x) {
    if( find(x) == end() ) return map_t::operator[](x) = x;
    Item y = x;
    while( map_t::operator[](y) != y ) y = map_t::operator[](y);
    return map_t::operator[](x) = y;
}

template<class Item, class Cmp>
void equivalence<Item,Cmp>::addrel(const Item& x,const Item& y) {
    Item z = operator[](x);
    Item t = operator[](y);

    if( cmp(z,t) ) map_t::operator[](z) = t;
              else map_t::operator[](t) = z;
}

template<class Item, class Cmp>
bool equivalence<Item,Cmp>::isequiv(const Item& x,const Item& y) {
    return operator[](x)==operator[](y);
}

template<class Item, class Cmp>
void equivalence<Item,Cmp>::update() {
    for(typename map_t::const_iterator mi=begin();mi!=end();++mi)
	operator[](mi->first);
}

template<class Item, class Cmp>
int equivalence<Item,Cmp>::classes() {
    int c = 0;
    for(typename map_t::const_iterator mi=begin();mi!=end();++mi)
	if(operator[](mi->first)==mi->first) ++c;
    return c;
}

template<class Item, class Cmp>
template<class eclass_t>
void equivalence<Item,Cmp>::get_eclasses(map<Item,eclass_t,Cmp>& C) {
    C.clear();
    for(typename map_t::const_iterator mi=begin();mi!=end();++mi)
        C[operator[](mi->first)].insert(mi->first);
}

/*
template<class Item, class Cmp>
template<class eclass_t>
void equivalence<Item,Cmp>::get_eclasses(map<Item,eclass_t,Cmp>& C) const {
    C.clear();
    for(typename map_t::const_iterator mi=begin();mi!=end();++mi)
        C[mi->second].insert(mi->first);
}
*/

#endif
