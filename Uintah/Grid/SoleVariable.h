
#ifndef UINTAH_HOMEBREW_SoleVariable_H
#define UINTAH_HOMEBREW_SoleVariable_H

#include "DataItem.h"
#include "TypeMismatchException.h"
#include <iostream> // TEMPORARY

template<class T>
class SoleVariable : public DataItem {
public:
    inline SoleVariable() {}
    inline SoleVariable(T value) : value(value) {}
    inline SoleVariable(const SoleVariable<T>& copy) : value(copy.value) {}
    virtual ~SoleVariable();

    static const TypeDescription* getTypeDescription();

    inline operator T () const {
	return value;
    }
    virtual void get(DataItem&) const;
    virtual SoleVariable<T>* clone() const;
    virtual void allocate(const Region*);
private:
    SoleVariable<T>& operator=(const SoleVariable<T>& copy);
    T value;
};

template<class T>
const TypeDescription* SoleVariable<T>::getTypeDescription()
{
    //cerr << "SoleVariable::getTypeDescription not done\n";
    return 0;
}

template<class T>
void SoleVariable<T>::get(DataItem& copy) const
{
    SoleVariable<T>* ref = dynamic_cast<SoleVariable<T>*>(&copy);
    if(!ref)
	throw TypeMismatchException("SoleVariable<T>");
    *ref = *this;
}

template<class T>
SoleVariable<T>::~SoleVariable()
{
}

template<class T>
SoleVariable<T>* SoleVariable<T>::clone() const
{
    return new SoleVariable<T>(*this);
}

template<class T>
SoleVariable<T>& SoleVariable<T>::operator=(const SoleVariable<T>& copy)
{
    value = copy.value;
    return *this;
}

template<class T>
void SoleVariable<T>::allocate(const Region*)
{
    throw TypeMismatchException("SoleVariable shouldn't use allocate");
}

#endif
