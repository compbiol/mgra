#ifndef PROPERTY_H_
#define PROPERTY_H_

/*
 *  Simple instance iterator property 
 */
#define DECLARE_ITERATOR(type, member, name, funcname) type name() {return this->member.funcname();}
#define DECLARE_CONST_ITERATOR(type, member, name, funcname) type name() const {return this->member.funcname();}

/*
 * Simple instance variable property
 */

#define DECLARE_GETTER(type, member, name) type get##name() const {return this->member;}
#define DECLARE_SETTER(type, member, name) void set##name(type v){this->member = v;}

#define DECLARE_PROPERTY(type,member,name) \
  DECLARE_GETTER(type,member,name)\
  DECLARE_SETTER(type,member,name)

/*
 * Property of member object
 */

#define DECLARE_DELEGATE_GETTER(type, member, name) type get##name() const {return this->member.get##name();}
#define DECLARE_DELEGATE_SETTER(type, member, name) void set##name(type v){this->member.set##name(v);}

#endif 
