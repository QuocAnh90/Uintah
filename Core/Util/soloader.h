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

// soloader.h written by Chris Moulding 11/98
// these functions are used to abstract the interface for 
// accessing shared libraries (.so for unix and .dll for windows)

#include <Core/share/share.h>
#include <string>

#ifdef _WIN32
#include <afxwin.h> // for LoadLibrary(), GetProcAddress() and HINSTANCE
typedef HINSTANCE LIBRARY_HANDLE;
#else
#include <dlfcn.h>   // for dlopen() and dlsym()
typedef void* LIBRARY_HANDLE;
#endif

/////////////////////////////////////////////////////////////////
// GetLibrarySymbolAddress()
// returns a pointer to the data or function called "symbolname"
// from within the shared library called "libname"

SCICORESHARE void* GetLibrarySymbolAddress(const char* libname, const char* symbolname);


/////////////////////////////////////////////////////////////////
// GetLibraryHandle()
// opens, and returns the handle to, the library module
// called "libname"

SCICORESHARE LIBRARY_HANDLE GetLibraryHandle(const char* libname);


/////////////////////////////////////////////////////////////////
// GetHandleSymbolAddress()
// returns a pointer to the data or function called "symbolname"
// from within the shared library with handle "handle"

SCICORESHARE void* GetHandleSymbolAddress(LIBRARY_HANDLE handle, const char* symbolname);


/////////////////////////////////////////////////////////////////
// CloseLibrary()
//
// disassociates the given library handle from the calling process.

SCICORESHARE void CloseLibrary(LIBRARY_HANDLE);


////////////////////////////////////////////////////////////////
// SOError()
// returns the last error generated by one of the above functions.

const char* SOError(bool =false);

LIBRARY_HANDLE FindLibInPath(const std::string& lib, const std::string& path);

