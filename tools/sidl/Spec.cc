
#include "Spec.h"
#include "SymbolTable.h"
#include <map>
#include <pair.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>

using std::cerr;
using std::find;
using std::map;
using std::string;
using std::vector;

Definition::Definition(const string& curfile, int lineno, const string& name)
    : curfile(curfile), lineno(lineno), name(name)
{
    symbols=0;
    emitted_declaration=false;
}

Definition::~Definition()
{
}

string Definition::fullname() const
{
    if(!symbols){
	cerr << "ERROR: Definition symboltable not set!\n";
	exit(1);
    }
    return symbols->fullname();
}

string CI::cppfullname(SymbolTable* forpackage) const
{
    if(!symbols){
	cerr << "ERROR: Definition symboltable not set!\n";
	exit(1);
    }
    if(forpackage == symbols->getParent())
	return name;
    else
	return symbols->cppfullname();
}

string CI::cppclassname() const
{
    return name;
}

SymbolTable* Definition::getSymbolTable() const
{
    return symbols;
}

Argument::Argument(Argument::Mode mode, Type* type)
    : mode(mode), type(type)
{
}

Argument::Argument(Argument::Mode mode, Type* type, const string& id)
    : mode(mode), type(type), id(id)
{
}

Argument::~Argument()
{
}

string Argument::fullsignature() const
{
    string s;
    switch(mode){
    case In:
	s="in ";
	break;
    case Out:
	s="out ";
	break;
    case InOut:
	s="inout ";
	break;
    }
    s+=type->fullname();
    if(id != "")
	s+=" "+id;
    return s;
}

Argument::Mode Argument::getMode() const
{
    return mode;
}

ArgumentList::ArgumentList()
{
}

ArgumentList::~ArgumentList()
{
    for(vector<Argument*>::iterator iter=list.begin(); iter != list.end();iter++){
	delete *iter;
    }
}

void ArgumentList::prepend(Argument* arg)
{
    list.insert(list.begin(), arg);
}

void ArgumentList::add(Argument* arg)
{
    list.push_back(arg);
}

std::vector<Argument*>& ArgumentList::getList()
{
    return list;
}

string ArgumentList::fullsignature() const
{
    std::string s="";
    int c=0;
    for(vector<Argument*>::const_iterator iter=list.begin(); iter != list.end(); iter++){
	if(c++ > 0)
	    s+=", ";
	s+=(*iter)->fullsignature();
    }
    return s;
}

CI::CI(const string& curfile, int lineno, const string& name,
       MethodList* methods)
    : Definition(curfile, lineno, name), mymethods(methods)
{
    vtable_base=-1234;
    parentclass=0;
}

CI::~CI()
{
    if(mymethods)
	delete mymethods;
}

Class::Class(const string& curfile, int lineno, const string& name,
	     ScopedName* class_extends, ScopedNameList* class_implements,
	     MethodList* methods)
    : CI(curfile, lineno, name, methods),
      class_extends(class_extends), class_implements(class_implements)
{
}

Class::Class(const string& curfile, int lineno, const string& name)
    : CI(curfile, lineno, name, 0),
      class_extends(0), class_implements(0)
{
}

Class::~Class()
{
    if(class_extends)
	delete class_extends;
    if(class_implements)
	delete class_implements;
}

void CI::gatherMethods(vector<Method*>& allmethods) const
{
    if(parentclass)
	parentclass->gatherMethods(allmethods);
    for(vector<Interface*>::const_iterator iter=parent_ifaces.begin();
	iter != parent_ifaces.end(); iter++){
	(*iter)->gatherMethods(allmethods);
    }
    mymethods->gatherMethods(allmethods);
}

void CI::gatherVtable(vector<Method*>& uniquemethods) const
{
    if(parentclass)
	parentclass->gatherVtable(uniquemethods);
    for(vector<Interface*>::const_iterator iter=parent_ifaces.begin();
	iter != parent_ifaces.end(); iter++){
	(*iter)->gatherVtable(uniquemethods);
    }
    mymethods->gatherVtable(uniquemethods);
}

vector<Method*>& CI::myMethods()
{
    return mymethods->getList();
}

void CI::gatherParents(std::vector<CI*>& folks) const
{
    if(parentclass)
	parentclass->gatherParents(folks);
    for(vector<Interface*>::const_iterator iter=parent_ifaces.begin();
	iter != parent_ifaces.end(); iter++){
	(*iter)->gatherParents(folks);
    }
    if(find(folks.begin(), folks.end(), this) == folks.end())
	folks.push_back(const_cast<CI*>(this));
}

void CI::gatherParentInterfaces(std::vector<Interface*>& folks) const
{
    if(parentclass)
	parentclass->gatherParentInterfaces(folks);
    for(vector<Interface*>::const_iterator iter=parent_ifaces.begin();
	iter != parent_ifaces.end(); iter++){
	(*iter)->gatherParentInterfaces(folks);
	if(find(folks.begin(), folks.end(), *iter) == folks.end())
	    folks.push_back(const_cast<Interface*>(*iter));
    }
}

MethodList* CI::getMethods() const
{
    return mymethods;
}

Class* Class::getParentClass() const
{
    return parentclass;
}

Class* Class::findParent(const string& str)
{
    Class* p=this;
    while(p){
	if(p->fullname() == str)
	    return p;
	p=p->parentclass;
    }
    return 0;
}

DefinitionList::DefinitionList()
{
}

DefinitionList::~DefinitionList()
{
    for(vector<Definition*>::iterator iter=list.begin();iter != list.end();iter++){
	delete *iter;
    }
}

void DefinitionList::add(Definition* def)
{
    list.push_back(def);
}

Interface::Interface(const string& curfile, int lineno, const string& id,
		     ScopedNameList* interface_extends, MethodList* methods)
    : CI(curfile, lineno, id, methods),
      interface_extends(interface_extends)
{
}

Interface::Interface(const string& curfile, int lineno, const string& id)
    : CI(curfile, lineno, id, 0), interface_extends(0)
{
}

Interface::~Interface()
{
    if(interface_extends)
	delete interface_extends;
}

Method::Method(const string& curfile, int lineno, Type* return_type,
	       const string& name, ArgumentList* args,
	       ScopedNameList* throws_clause)
    : curfile(curfile), lineno(lineno), return_type(return_type),
      name(name), args(args),
      throws_clause(throws_clause), modifier(Unknown)
{
    myclass=0;
    myinterface=0;
    checked=false;
    handlerNum=-1;
}

Method::~Method()
{
    delete args;
    delete throws_clause;
}

void Method::setModifier(Method::Modifier newmodifier)
{
    modifier=newmodifier;
}

Method::Modifier Method::getModifier() const
{
    return modifier;
}

void Method::setClass(Class *c)
{
    myclass=c;
}

void Method::setInterface(Interface *c)
{
    myinterface=c;
}

bool Method::getChecked() const
{
    return checked;
}

const string& Method::getCurfile() const
{
    return curfile;
}

int Method::getLineno() const
{
    return lineno;
}

string Method::fullname() const
{
    string n;
    if(myclass)
	n=myclass->fullname();
    else if(myinterface)
	n=myinterface->fullname();
    else {
	cerr << "ERROR: Method does not have an associated class or interface\n";
	exit(1);
    }
    return n+"."+name;
}

string Method::fullsignature() const
{
    string s;
    switch(modifier){
    case Abstract:
	s="abstract ";
	break;
    case Final:
	s="final ";
	break;
    case Static:
	s="static ";
	break;
    case Unknown:
	s="";
	break;
    }
    s+=return_type->fullname()+" "+fullname()+"("+args->fullsignature()+")";
    if(throws_clause){
	s+="throws "+throws_clause->fullsignature();
    }
    return s;
}

MethodList::MethodList()
{
}

MethodList::~MethodList()
{
    for(vector<Method*>::iterator iter=list.begin(); iter != list.end(); iter++){
	delete *iter;
    }
}

void MethodList::add(Method* method)
{
    list.push_back(method);
}

void MethodList::setClass(Class* c)
{
    for(vector<Method*>::iterator iter=list.begin();iter != list.end();iter++){
	(*iter)->setClass(c);
    }
}

void MethodList::setInterface(Interface* c)
{
    for(vector<Method*>::iterator iter=list.begin();iter != list.end();iter++){
	(*iter)->setInterface(c);
    }
}

vector<Method*>& MethodList::getList()
{
    return list;
}

void MethodList::gatherMethods(vector<Method*>& allmethods) const
{
    for(vector<Method*>::const_iterator iter=list.begin();iter != list.end();iter++){
	allmethods.push_back(*iter);
    }
}

void MethodList::gatherVtable(vector<Method*>& uniquemethods) const
{
    // Yuck - O(N^2)
    for(vector<Method*>::const_iterator iter1=list.begin();iter1 != list.end();iter1++){
	vector<Method*>::iterator iter2=uniquemethods.begin();;
	for(;iter2 != uniquemethods.end(); iter2++){
	    if((*iter1)->matches(*iter2, Method::TypesOnly))
		break;
	}
	if(iter2 == uniquemethods.end())
	    uniquemethods.push_back(*iter1);
    }
}

Package::Package(const string& fileno, int lineno, const string& name,
		 DefinitionList* definition)
    : Definition(fileno, lineno, name), name(name), definition(definition)
{
}

Package::~Package()
{
    delete definition;
}

ScopedName::ScopedName()
    : leading_dot(false)
{
    sym=0;
}

ScopedName::ScopedName(const string& s1, const string& s2)
{
    leading_dot=true;
    add(s1);
    add(s2);
    sym=0;
}

ScopedName::~ScopedName()
{
}

void ScopedName::prepend(const string& name)
{
    names.insert(names.begin(), name);
}

void ScopedName::add(const string& name)
{
    names.push_back(name);
}

void ScopedName::set_leading_dot(bool dot)
{
    leading_dot=dot;
}

string ScopedName::getName() const
{
    string n;
    if(leading_dot)
	n+=".";
    bool first=true;
    for(vector<string>::const_iterator iter=names.begin();iter != names.end();iter++){
	if(!first)
	    n+=".";
	n+=*iter;
	first=false;
    }
    return n;
}

Symbol* ScopedName::getSymbol() const
{
    return sym;
}

bool ScopedName::getLeadingDot() const
{
    return leading_dot;
}

const string& ScopedName::name(int idx) const
{
    return names[idx];
}

int ScopedName::nnames() const
{
    return names.size();
}

string ScopedName::fullname() const
{
    if(!sym){
	cerr << "ERROR: Symbol not bound: " << getName() << " - " << this << "\n";
	exit(1);
    }
    return sym->fullname();
}

string ScopedName::cppfullname(SymbolTable* forstab) const
{
    if(!sym){
	cerr << "ERROR: Symbol not bound: " << getName() << " - " << this << "\n";
	exit(1);
    }
    return sym->cppfullname(forstab);
}

ScopedNameList::ScopedNameList()
{
}

ScopedNameList::~ScopedNameList()
{
    for(vector<ScopedName*>::iterator iter=list.begin(); iter != list.end(); iter++){
	delete *iter;
    }
}

void ScopedNameList::prepend(ScopedName* name)
{
    list.insert(list.begin(), name);
}

void ScopedNameList::add(ScopedName* name)
{
    list.push_back(name);
}

vector<ScopedName*> const& ScopedNameList::getList() const
{
    return list;
}

std::string ScopedNameList::fullsignature() const
{
    string s="";
    int c=0;
    for(vector<ScopedName*>::const_iterator iter=list.begin(); iter != list.end(); iter++){
	if(c++ > 0)
	    s+=", ";
	s+=(*iter)->fullname();
    }
    return s;
}

Specification::Specification()
{
}

Specification::~Specification()
{
    for(vector<DefinitionList*>::iterator iter=list.begin();iter != list.end(); iter++){
	delete *iter;
    }
}

void Specification::add(DefinitionList* def)
{
    list.push_back(def);
}

Type* Type::voidtype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("void", "error");
    return t;
}

Type* Type::booltype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("bool", "byte");
    return t;
}

Type* Type::chartype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("char", "char");
    return t;
}

Type* Type::dcomplextype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("dcomplex", "special");
    return t;
}

Type* Type::doubletype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("double", "double");
    return t;
}

Type* Type::fcomplextype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("fcomplex", "special");
    return t;
}

Type* Type::floattype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("float", "float");
    return t;
}

Type* Type::inttype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("int", "int");
    return t;
}

Type* Type::stringtype()
{
    static Type* t;
    if(!t)
	t=new BuiltinType("string", "special");
    return t;
}

static map<pair<Type*, int>, Type*> arrays;

Type* Type::arraytype(Type* t, int dim)
{
    map<pair<Type*, int>, Type*>::iterator iter=arrays.find(pair<Type*, int>(t, dim));
    if(iter == arrays.end()){
	return (arrays[pair<Type*, int>(t, dim)]=new ArrayType(t, dim));
    } else {
	return iter->second;
    }
}

Type::Type()
{
}

Type::~Type()
{
}

BuiltinType::BuiltinType(const string& cname, const string& nexusname)
    : cname(cname), nexusname(nexusname)
{
}

BuiltinType::~BuiltinType()
{
}

std::string BuiltinType::fullname() const
{
    return cname;
}

bool BuiltinType::isvoid() const
{
    if(cname == "void")
	return true;
    else
	return false;
}

NamedType::NamedType(const string& curfile, int lineno, ScopedName* name)
    : curfile(curfile), lineno(lineno), name(name)
{
}

NamedType::~NamedType()
{
}

std::string NamedType::fullname() const
{
    return name->fullname();
}

bool NamedType::isvoid() const
{
    return false;
}

ArrayType::ArrayType(Type* t, int dim)
    : subtype(t), dim(dim)
{
}

ArrayType::~ArrayType()
{
}

std::string ArrayType::fullname() const
{
    std::ostringstream o;
    o << "array<" << subtype->fullname() << ", " << dim << ">";
    return o.str();
}

bool ArrayType::isvoid() const
{
    return false;
}


//
// $Log$
// Revision 1.5  1999/09/21 06:13:01  sparker
// Fixed bugs in multiple inheritance
// Added round-trip optimization
// To support this, we store Startpoint* in the endpoint instead of the
//    object final type.
//
// Revision 1.4  1999/09/17 05:07:25  sparker
// Added nexus code generation capability
//
// Revision 1.3  1999/08/30 20:19:29  sparker
// Updates to compile with -LANG:std on SGI
// Other linux/irix porting oscillations
//
// Revision 1.2  1999/08/30 17:39:54  sparker
// Updates to configure script:
//  rebuild configure if configure.in changes (Bug #35)
//  Fixed rule for rebuilding Makefile from Makefile.in (Bug #36)
//  Rerun configure if configure changes (Bug #37)
//  Don't build Makefiles for modules that aren't --enabled (Bug #49)
//  Updated Makfiles to build sidl and Component if --enable-parallel
// Updates to sidl code to compile on linux
// Imported PIDL code
// Created top-level Component directory
// Added ProcessManager class - a simpler interface to fork/exec (not finished)
//
//
