#ifndef PROPERTY_H_
#define PROPERTY_H_

/*
 *  Simple instance iterator property 
 */
#define DECLARE_ITERATOR(type, member, name, funcname) type name() {return this->member.funcname();}
#define DECLARE_CONST_ITERATOR(type, member, name, funcname) type name() const {return this->member.funcname();}

/*
 * Delegators
 */
#define DECLARE_DELEGATE_VOID_METHOD(member, name, funcname) void name() {this->member.funcname();}

#define DECLARE_DELEGATE_METHOD(type, member, name, funcname) type name() {return this->member.funcname();}
#define DECLARE_DELEGATE_CONST_METHOD(type, member, name, funcname) type name() const {return this->member.funcname();}

#define DECLARE_DELEGATE_PARAM_METHOD(type, member, name, typep, funcname) type name(typep v) {return this->member.funcname(v);}
#define DECLARE_DELEGATE_PARAM_CONST_METHOD(type, member, name, typep, funcname) type name(typep v) const {return this->member.funcname(v);}


//#define DECLARE_DELEGATE_VOID_METHOD(member, name, funcname) void name() {this->member.funcname();}
//#define DECLARE_DELEGATE_VOID_CONST_METHOD(member, name, funcname) void name() const {this->member.funcname();}


/*
 * Simple instance variable property
 */
#define DECLARE_GETTER(type, member, name) type get_##name() const {return this->member;}
#define DECLARE_SETTER(type, member, name) void set_##name(type v){this->member = v;}

//#define DECLARE_DELEGATE_GETTER(type, member, name) type get_##name() const {return this->member;}


/*
 *
 */
#define DECLARE_DEFINED(type, check, name) bool name(type const & v) const {return check;}

/*
 * While not usedget_partgraph
 */


/*
#define DECLARE_PROPERTY(type,member,name) \
  DECLARE_GETTER(type,member,name)\
  DECLARE_SETTER(type,member,name)
*/
/*
 * Property of member object
 */

//#define DECLARE_DELEGATE_GETTER(type, member, name) type get##name() const {return this->member.get##name();}
//#define DECLARE_DELEGATE_SETTER(type, member, name) void set##name(type v){this->member.set##name(v);}


//#define DECLARE_DELEGATE_METHOD(type, member, name, funcname) type name() {return this->member.funcname();}
//#define DECLARE_DELEGATE_VOID_METHOD(type, member, name) void set##name(type v){this->member = v;}

#endif 
