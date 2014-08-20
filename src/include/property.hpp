#ifndef PROPERTY_HPP
#define PROPERTY_HPP

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


/*
 * Simple instance variable property
 */
#define DECLARE_GETTER(type, member, name) type get_##name() const {return this->member;}
#define DECLARE_SETTER(type, member, name) void set_##name(type v){this->member = v;}

#endif 
