
/*
 *  String.cc: implementation of String utility class
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   Feb. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */


#include <Classlib/String.h>
#include <Classlib/Assert.h>
#include <Classlib/Persistent.h>
#include <Classlib/TrivialAllocator.h>
#include <Malloc/Allocator.h>
#include <iostream.h>
#include <stdio.h>
#include <string.h>

#ifdef BROKEN
static TrivialAllocator* srep_alloc=0;

inline void* clString::srep::operator new(size_t)
{
    if(!srep_alloc){
	srep_alloc=scinew TrivialAllocator(sizeof(clString::srep));
	lock=scinew Mutex;
    }
    lock->lock();
    void* p=srep_alloc->alloc();
    lock->unlock();
    return p;
}

inline void clString::srep::operator delete(void* rp, size_t)
{
    lock->lock();
    srep_alloc->free(rp);
    lock->unlock();
}
#endif

#define inline
#include <Classlib/String.icc>
#undef inline

clString::clString(const char* s)
{
    p=scinew srep;
    int len=strlen(s);
    p->s=scinew char[len+1];
    strcpy(p->s,s);
}

clString clString::operator+(const clString& str) const
{
    int newlen=(p?strlen(p->s):0) + (str.p?strlen(str.p->s):0);
    char* ns=scinew char[newlen+1];
    if(p && p->s)strcpy(ns, p->s);
    else ns[0]=0;
    if(str.p && str.p->s)strcat(ns, str.p->s);
    return clString(0, ns);
}

clString clString::operator+(const char* c) const
{
    int newlen=(p?strlen(p->s):0)+strlen(c);
    char* ns=scinew char[newlen+1];
    if(p && p->s)strcpy(ns, p->s);
    else ns[0]=0;
    strcat(ns, c);
    return clString(0, ns);
}

clString operator+(const char* c, const clString& str)
{
    int newlen=(str.p?strlen(str.p->s):0)+strlen(c);
    char* ns=scinew char[newlen+1];
    strcpy(ns, c);
    if(str.p && str.p->s)strcat(ns, str.p->s);
    return clString(0, ns);
}

int clString::get_double(double& x) const
{
    return (p && p->s)?sscanf(p->s, "%lf", &x)==1:0;
}

int clString::get_int(int& x) const
{
    return (p && p->s)?sscanf(p->s, "%d", &x)==1:0;
}

clString to_string(int n)
{
    char s[50];
    sprintf(s,"%d",n);
    return clString(s);
}

clString to_string(double d)
{
    char s[50];
    sprintf(s,"%g",d);
    return clString(s);
}

ostream& operator<<(ostream& s, const clString& str)
{
    return s << ((str.p && str.p->s)?str.p->s:"");
}

istream& operator>>(istream& s, clString& str)
{
    char* buf=scinew char[1000];
    s.get(buf,1000,'\n');
#ifdef broken
    char c;
    if(cin.get(c) && c!='\n'){
	// Longer than 1000...
	int grow=1;
	int size=1000;
	while(grow){
	    int newsize=size << 1; /* Double size... */
	    char* p=scinew char[newsize];
	    strncpy(p, buf, size);
	    s.get(buf+size,size,'\n');
	    if(cin.get(c) && c!='\n'){
		grow=1;
	    } else {
		grow=0;
	    }
	    delete[] buf;
	    buf=p;
	    size=newsize;
	}
    }
#endif
    str=buf; // Uses operator=
    delete[] buf;
    return s;
}

int clString::index(const char match) const
{
    if(!p || !p->s)return -1;
    int i=0;
    char* pp=p->s;
    while(*pp){
	if(*pp == match)return i;
	i++;
	pp++;
    }
    return -1;
}

clString clString::substr(int start, int length)
{
    ASSERT(p != 0);
    int len=strlen(p->s);
    ASSERTRANGE(start, 0, len);
    int l=length==-1?len-start:length;
    ASSERTRANGE(start+l, 0, len+1);
    char* tmp=scinew char[l+1];
    int i;
    for(i=0;i<l;i++){
	tmp[i]=p->s[i+start];
    }
    tmp[i]='\0';
    clString rstr(0, tmp);
    return rstr;
}

int clString::hash(int hash_size) const
{
    if(!p || !p->s)return 0;
    char* pp=p->s;
    int sum=0;
    while(*pp){
	sum=(sum << 2) ^ (sum >> 2) ^ *pp;
	pp++;
    }
    sum=sum<0?-sum:sum;
    return sum%hash_size;
}

clString basename(const clString& str)
{
    ASSERT(str.p && str.p->s);
    char* pp=str.p->s;
    char* last_slash=pp;
    while(*pp){
	if(*pp=='/')last_slash=pp;
	pp++;
    }
    return clString(last_slash+1);
}

clString& clString::operator+=(char c)
{
    if(p){
	if(p->n != 1){
	    // detach...
	    srep* oldp=p;
	    p=scinew srep;
	    int len=strlen(oldp->s);
	    p->s=scinew char[len+2];
	    strcpy(p->s, oldp->s);
	    p->s[len]=c;
	    p->s[len+1]=0;
	    oldp->n--;
	    p->n=1;
	} else {
	    char* olds=p->s;
	    int len=strlen(olds);
	    p->s=scinew char[len+2];
	    strcpy(p->s, olds);
	    p->s[len]=c;
	    p->s[len+1]=0;
	    delete[] olds;
	}
    } else {
	p=scinew srep;
	p->n=1;
	p->s=scinew char[2];
	p->s[0]=c;
	p->s[1]=0;
    }
    return *this;
}

clString& clString::operator+=(const clString& str)
{
    int newlen=(p?strlen(p->s):0)+(str.p?strlen(str.p->s):0);
    char* ns=scinew char[newlen+1];
    if(p && p->s)strcpy(ns, p->s);
    else ns[0]=0;
    if(str.p && str.p->s)strcat(ns, str.p->s);
    if(p && p->n > 1){
	if(p)p->n--;
	p=scinew srep;
    } else {
	if(p && p->s)delete[] p->s;
	if(!p)
	    p=scinew srep;
    }
    p->s=ns;
    return *this;
}

#include <Tester/RigorousTest.h>

void clString::test_rigorous(RigorousTest* __test)
{
    TEST(1==0);
    
    clString s0("hi");
    clString s1(s0+" ");
    clString s2("there");
    TEST(s0<s2);
    TEST(s2>s0);
    clString s3(s1+s2);
    TEST(s3=="hi there");
    clString s4(s3);
    TEST(s1 != s4);
    clString s5(s4());
    TEST(s5==s4);
    TEST(s5.len()==8);
    TEST(s3(0)=='h');
    TEST(s3(2)==' ');
    TEST(s3(7)=='e');
    s3="";
    int i;
    int n=8;
    for(i=0;i<n;i++){
	s3=s3+"0123456789";
    }
    for(i=0;i<n;i++){
	s3=s3+"0123456789";
    }
    n=n+n;
    for(i=0;i<3;i++){
	s3=s3+s3;
    }
    n<<=3;
    for(i=0;i<3;i++){
	s3+=s3;
    }
    n<<=3;
    TEST(s3.len()==n*10);
    int idx=0;
    for(i=0;i<n;i++){
	for(int j=0;j<10;j++){
	    TEST(s3(idx++)==j+'0');
	}
    }
}

#include <Tester/PerfTest.h>

/*
 * Performance tests
 */
void clString::test_performance(PerfTest* __pt) {
    PERFTEST("baseline") {
    }

    /*
     * Test CTORs
     */
    PERFTEST("CTOR/DTOR - empty") {
	clString str;
    }
    PERFTEST("CTOR/DTOR - with small string") {
	clString str("asdf");
    }
    static char* large_string="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    PERFTEST("CTOR/DTOR - with large string") {
	clString str(large_string);
    }
    clString s1("asdf");
    PERFTEST("CTOR - copy small string") {
	clString str(s1);
    }
    clString s2(large_string);
    PERFTEST("CTOR - copy large string") {
	clString str(s2);
    }
    
    /*
     * Various operators.
     */
#if 0
    clString s3;
    PERFTEST("operator= - small string from clString") {
	s3=s1;
    }
    clString s4;
    PERFTEST("operator= - large string from clString") {
	s4=s2;
    }
    PERFTEST("operator= - small string from char*") {
	s3="asdf";
    }
    PERFTEST("operator= - large string from char*") {
	s4=large_string;
    }
    PERFTEST("operator+ - pairwise") {
	s4=s1+s2;
    }
#endif
    PERFTEST("compound + - big") {
	clString str(s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2);
    }
    PERFTEST("compound + - flatten big") {
	clString str(s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2+s2);
	str();
    }
}
