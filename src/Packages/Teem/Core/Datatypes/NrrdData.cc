/*
  The contents of this file are subject to the University of Utah Public
  License (the "License"); you may not use this file except in compliance
  with the License.
  
  Software distributed under the License is distributed on an "AS IS"
  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
  License for the specific language governing rights and limitations under
  the License.
  
  The Original Source Code is SCIRun, released March 12, 2001.
  
  The Original Source Code was developed by the University of Utah.
  Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
  University of Utah. All Rights Reserved.
*/

// NrrdData.cc - Interface to Gordon's Nrrd package
//
//  Written by:
//   David Weinstein
//   School of Computing
//   University of Utah
//   February 2001
//
//  Copyright (C) 2001 SCI Institute

#include <Teem/Core/Datatypes/NrrdData.h>
#include <Core/Malloc/Allocator.h>
#include <iostream>

using std::cerr;

namespace SCITeem {

static Persistent* make_NrrdData() {
  return scinew NrrdData;
}

PersistentTypeID NrrdData::type_id("NrrdData", "Datatype", make_NrrdData);

NrrdData::NrrdData(bool owned) : 
  nrrd(nrrdNew()),
  data_owned_(owned)
{}

NrrdData::NrrdData(const NrrdData &copy) :
  fname(copy.fname) 
{
  nrrd = nrrdNew();
  nrrdCopy(nrrd, copy.nrrd);
}

NrrdData::~NrrdData() {
  if(data_owned_) {
    nrrdNuke(nrrd);
  } else {
    nrrdNix(nrrd);
  }
}

// This needs to parse axis 0 and see if the label is tuple as well...
bool
NrrdData::is_sci_nrrd() const 
{
  return (originating_field_.get_rep() != 0);
}

void 
NrrdData::copy_sci_data(const NrrdData &cp)
{
  originating_field_ = cp.originating_field_;
}

#define NRRDDATA_VERSION 1

//////////
// PIO for NrrdData objects
void NrrdData::io(Piostream& stream) {
  /*  int version = */ stream.begin_class("NrrdData", NRRDDATA_VERSION);

  if (stream.reading()) {
    Pio(stream, fname);
    if (!(nrrdLoad(nrrd=nrrdNew(), strdup(fname.c_str())))) {
      char *err = biffGet(NRRD);
      cerr << "Error reading nrrd "<<fname<<": "<<err<<"\n";
      free(err);
      biffDone(NRRD);
      return;
    }
    fname="";
  } else { // writing
    if (fname == "") {   // if fname wasn't set up stream, just append .nrrd
      fname = stream.file_name + string(".nrrd");
    }
    Pio(stream, fname);
    if (nrrdSave(strdup(fname.c_str()), nrrd, 0)) {
      char *err = biffGet(NRRD);      
      cerr << "Error writing nrrd "<<fname<<": "<<err<<"\n";
      free(err);
      biffDone(NRRD);
      return;
    }
  }
  stream.end_class();
}


template <>
unsigned int get_nrrd_type<char>() {
  return nrrdTypeChar;
}


template <>
unsigned int get_nrrd_type<unsigned char>()
{
  return nrrdTypeUChar;
}

template <>
unsigned int get_nrrd_type<short>()
{
  return nrrdTypeShort;
}

template <>
unsigned int get_nrrd_type<unsigned short>()
{
  return nrrdTypeUShort;
}

template <>
unsigned int get_nrrd_type<int>()
{
  return nrrdTypeInt;
}

template <>
unsigned int get_nrrd_type<unsigned int>()
{
  return nrrdTypeUInt;
}

template <>
unsigned int get_nrrd_type<long long>()
{
  return nrrdTypeLLong;
}

template <>
unsigned int get_nrrd_type<unsigned long long>()
{
  return nrrdTypeULLong;
}

template <>
unsigned int get_nrrd_type<float>()
{
  return nrrdTypeFloat;
}



}  // end namespace SCITeem
