#ifndef UINTAH_HOMEBREW_PARTICLEVARIABLE_H
#define UINTAH_HOMEBREW_PARTICLEVARIABLE_H

#include <Core/Util/FancyAssert.h>
#include <Core/Exceptions/ErrnoException.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Util/Assert.h>
#include <Core/Malloc/Allocator.h>
#include <Dataflow/XMLUtil/XMLUtil.h>
#include <Packages/Uintah/Core/Grid/ParticleVariableBase.h>
#include <Packages/Uintah/CCA/Ports/InputContext.h>
#include <Packages/Uintah/CCA/Ports/OutputContext.h>
#include <Packages/Uintah/Core/Exceptions/TypeMismatchException.h>
#include <Packages/Uintah/Core/Grid/ParticleData.h>
#include <Packages/Uintah/Core/Grid/ParticleSubset.h>
#include <Packages/Uintah/Core/Grid/TypeDescription.h>
#include <Packages/Uintah/Core/Grid/TypeUtils.h>
#include <Packages/Uintah/Core/Parallel/ProcessorGroup.h>

#include <unistd.h>
#include <errno.h>

namespace Uintah {

using namespace SCIRun;

class TypeDescription;

/**************************************

CLASS
   ParticleVariable
   
   Short description...

GENERAL INFORMATION

   ParticleVariable.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   Copyright (C) 2000 SCI Group

KEYWORDS
   Particle_Variable

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

   template<class T> class ParticleVariable : public ParticleVariableBase {
   public:
      ParticleVariable();
      virtual ~ParticleVariable();
      ParticleVariable(ParticleSubset* pset);
      ParticleVariable(ParticleData<T>*, ParticleSubset* pset);
      ParticleVariable(const ParticleVariable<T>&);
      
      ParticleVariable<T>& operator=(const ParticleVariable<T>&);

      //////////
      // Insert Documentation Here:
      static const TypeDescription* getTypeDescription();

      //////////
      // Insert Documentation Here:
      void resync() {
	 d_pdata->resize(getParticleSet()->numParticles());
      }
      
      //////////
      // Insert Documentation Here:
      virtual ParticleVariableBase* clone() const;
      virtual ParticleVariableBase* cloneSubset(ParticleSubset*) const;
      
      //////////
      // Insert Documentation Here:
      inline T& operator[](particleIndex idx) {
	 ASSERTRANGE(idx, 0, (particleIndex)d_pdata->data.size());
	 return d_pdata->data[idx];
      }
      
      //////////
      // Insert Documentation Here:
      inline const T& operator[](particleIndex idx) const {
	 ASSERTRANGE(idx, 0, (particleIndex)d_pdata->data.size());
	 return d_pdata->data[idx];
      }
      
      virtual void copyPointer(const ParticleVariableBase&);
      virtual void allocate(ParticleSubset*);
      virtual void gather(ParticleSubset* dest,
			  std::vector<ParticleSubset*> subsets,
			  std::vector<ParticleVariableBase*> srcs,
			  particleIndex extra = 0);
      virtual void unpackMPI(void* buf, int bufsize, int* bufpos,
			     const ProcessorGroup* pg, int start, int n);
      virtual void packMPI(void* buf, int bufsize, int* bufpos,
			   const ProcessorGroup* pg, int start, int n);
      virtual void packsizeMPI(int* bufpos,
			       const ProcessorGroup* pg, int start, int n);
      virtual void emit(OutputContext&);
      virtual void read(InputContext&);
      virtual void* getBasePointer();
      virtual const TypeDescription* virtualGetTypeDescription() const;
   private:
      
      //////////
      // Insert Documentation Here:
      ParticleData<T>* d_pdata;
      
      static TypeDescription::Register registerMe;
      static Variable* maker();
   };

   template<class T>
      TypeDescription::Register ParticleVariable<T>::registerMe(getTypeDescription());

   template<class T>
      const TypeDescription*
      ParticleVariable<T>::getTypeDescription()
      {
	 static TypeDescription* td;
	 if(!td){
	    td = scinew TypeDescription(TypeDescription::ParticleVariable,
					"ParticleVariable", &maker,
					fun_getTypeDescription((T*)0));
	 }
	 return td;
      }
   
   template<class T>
      Variable*
      ParticleVariable<T>::maker()
      {
	 return scinew ParticleVariable<T>();
      }
   
   template<class T>
      ParticleVariable<T>::ParticleVariable()
      : ParticleVariableBase(0), d_pdata(0)
      {
      }
   
   template<class T>
      ParticleVariable<T>::~ParticleVariable()
      {
	 if(d_pdata && d_pdata->removeReference())
	    delete d_pdata;
      }
   
   template<class T>
      ParticleVariable<T>::ParticleVariable(ParticleSubset* pset)
      : ParticleVariableBase(pset)
      {
	 d_pdata=scinew ParticleData<T>(pset->getParticleSet()->numParticles());
	 d_pdata->addReference();
      }
   
   template<class T>
      void ParticleVariable<T>::allocate(ParticleSubset* pset)
      {
	 if(d_pdata && d_pdata->removeReference())
	    delete d_pdata;
	 if(d_pset && d_pset->removeReference())
	    delete d_pset;

	 d_pset=pset;
	 d_pset->addReference();
	 d_pdata=scinew ParticleData<T>(pset->getParticleSet()->numParticles());
	 d_pdata->addReference();
      }
   
   template<class T>
      ParticleVariableBase*
      ParticleVariable<T>::clone() const
      {
	 return scinew ParticleVariable<T>(*this);
      }
   
   template<class T>
      ParticleVariableBase*
      ParticleVariable<T>::cloneSubset(ParticleSubset* pset) const
      {
	 return scinew ParticleVariable<T>(d_pdata, pset);
      }

   template<class T>
      ParticleVariable<T>::ParticleVariable(ParticleData<T>* pdata,
					    ParticleSubset* pset)
      : ParticleVariableBase(pset), d_pdata(pdata)
      {
	 if(d_pdata)
	    d_pdata->addReference();
      }
   
   template<class T>
      ParticleVariable<T>::ParticleVariable(const ParticleVariable<T>& copy)
      : ParticleVariableBase(copy), d_pdata(copy.d_pdata)
      {
	 if(d_pdata)
	    d_pdata->addReference();
      }
   
   template<class T>
      ParticleVariable<T>&
      ParticleVariable<T>::operator=(const ParticleVariable<T>& copy)
      {
	 if(this != &copy){
	    ParticleVariableBase::operator=(copy);
	    if(d_pdata && d_pdata->removeReference())
	       delete d_pdata;
	    d_pdata = copy.d_pdata;
	    if(d_pdata)
	       d_pdata->addReference();
	 }
	 return *this;
      }
   
   template<class T>
      void
      ParticleVariable<T>::copyPointer(const ParticleVariableBase& copy)
      {
	 const ParticleVariable<T>* c = dynamic_cast<const ParticleVariable<T>* >(&copy);
	 if(!c)
	    throw TypeMismatchException("Type mismatch in particle variable");
	 *this = *c;
      }

   template<class T>
      void
      ParticleVariable<T>::gather(ParticleSubset* pset,
				  std::vector<ParticleSubset*> subsets,
				  std::vector<ParticleVariableBase*> srcs,
				  particleIndex extra)
      {
	 if(d_pdata && d_pdata->removeReference())
	    delete d_pdata;
	 if(d_pset && d_pset->removeReference())
	    delete d_pset;
	 d_pset = pset;
	 pset->addReference();
	 d_pdata=scinew ParticleData<T>(pset->getParticleSet()->numParticles());
	 d_pdata->addReference();
	 ASSERTEQ(subsets.size(), srcs.size());
	 ParticleSubset::iterator dstiter = pset->begin();
	 for(int i=0;i<(int)subsets.size();i++){
	    ParticleVariable<T>* srcptr = dynamic_cast<ParticleVariable<T>*>(srcs[i]);
	    if(!srcptr)
	       throw TypeMismatchException("Type mismatch in ParticleVariable::gather");
	    ParticleVariable<T>& src = *srcptr;
	    ParticleSubset* subset = subsets[i];
	    for(ParticleSubset::iterator srciter = subset->begin();
		srciter != subset->end(); srciter++){
	       (*this)[*dstiter] = src[*srciter];
	       dstiter++;
	    }
	 }
	 ASSERT(dstiter+extra == pset->end());
      }

   template<class T>
      void
      ParticleVariable<T>::emit(OutputContext& oc)
      {
	 const TypeDescription* td = fun_getTypeDescription((T*)0);
	 appendElement(oc.varnode, "numParticles", d_pset->numParticles());
	 if(td->isFlat()){
	    // This could be optimized...
	    ParticleSubset::iterator iter = d_pset->begin();
	    while(iter != d_pset->end()){
	       particleIndex start = *iter;
	       iter++;
	       particleIndex end = start+1;
	       while(iter != d_pset->end() && *iter == end) {
		  end++;
		  iter++;
	       }
	       ssize_t size = (ssize_t)(sizeof(T)*(end-start));
	       ssize_t s=write(oc.fd, &(*this)[start], size);
	       if(size != s)
		  throw ErrnoException("ParticleVariable::emit (write call)", errno);
	       oc.cur+=size;
	    }
	 } else {
	    throw InternalError("Cannot yet write non-flat objects!\n");
	 }
      }

   template<class T>
      void
      ParticleVariable<T>::read(InputContext& ic)
      {
	 const TypeDescription* td = fun_getTypeDescription((T*)0);
	 if(td->isFlat()){
	    // This could be optimized...
	    ParticleSubset::iterator iter = d_pset->begin();
	    while(iter != d_pset->end()){
	       particleIndex start = *iter;
	       iter++;
	       particleIndex end = start+1;
	       while(iter != d_pset->end() && *iter == end) {
		  end++;
		  iter++;
	       }
	       ssize_t size = (ssize_t)(sizeof(T)*(end-start));
	       ssize_t s=::read(ic.fd, &(*this)[start], size);
	       if(size != s)
		  throw ErrnoException("ParticleVariable::emit (write call)", errno);
	       ic.cur+=size;
	    }
	 } else {
	    throw InternalError("Cannot yet write non-flat objects!\n");
	 }
      }

   template<class T>
      void*
      ParticleVariable<T>::getBasePointer()
      {
	 return &d_pdata->data[0];
      }

   template<class T>
      const TypeDescription*
      ParticleVariable<T>::virtualGetTypeDescription() const
      {
	 return getTypeDescription();
      }
   
   template<class T>
      void
      ParticleVariable<T>::unpackMPI(void* buf, int bufsize, int* bufpos,
				     const ProcessorGroup* pg,
				     int start, int n)
      {
	 // This should be fixed for variable sized types!
	 const TypeDescription* td = getTypeDescription()->getSubType();
	 if(td->isFlat()){
	    ParticleSubset::iterator beg = d_pset->seek(start);
	    ParticleSubset::iterator end = d_pset->seek(start+n);
	    for(ParticleSubset::iterator iter = beg; iter != end; iter++){
	       MPI_Unpack(buf, bufsize, bufpos,
			  &d_pdata->data[*iter], 1, td->getMPIType(),
			  pg->getComm());
	    }
	 } else {
	    throw InternalError("packMPI not finished\n");
	 }
      }

   template<class T>
      void
      ParticleVariable<T>::packMPI(void* buf, int bufsize, int* bufpos,
			   const ProcessorGroup* pg, particleIndex start, int n)
      {
	 // This should be fixed for variable sized types!
	 const TypeDescription* td = getTypeDescription()->getSubType();
	 if(td->isFlat()){
	    ParticleSubset::iterator beg = d_pset->seek(start);
	    ParticleSubset::iterator end = d_pset->seek(start+n);
	    for(ParticleSubset::iterator iter = beg; iter != end; iter++){
	       MPI_Pack(&d_pdata->data[*iter], 1, td->getMPIType(),
			buf, bufsize, bufpos, pg->getComm());
	    }
	 } else {
	    throw InternalError("packMPI not finished\n");
	 }
      }

   template<class T>
      void
      ParticleVariable<T>::packsizeMPI(int* bufpos,
			       const ProcessorGroup* pg, int, int n)
      {
	 // This should be fixed for variable sized types!
	 const TypeDescription* td = getTypeDescription()->getSubType();
	 if(td->isFlat()){
	    int size;
	    MPI_Pack_size(n, td->getMPIType(), pg->getComm(), &size);
	    (*bufpos)+= size;
	 } else {
	    throw InternalError("packsizeMPI not finished\n");
	 }
      }
} // End namespace Uintah



#endif
