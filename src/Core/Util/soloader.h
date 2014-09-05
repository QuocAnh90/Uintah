/*
 * The MIT License
 *
 * Copyright (c) 1997-2014 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

// soloader.h written by Chris Moulding 11/98
// these functions are used to abstract the interface for 
// accessing shared libraries (.so for unix and .dll for windows)

#ifndef SOLOADER_H_
#define SOLOADER_H_

#include <string>
#include <dlfcn.h>   // for dlopen() and dlsym()

typedef void* LIBRARY_HANDLE;

/////////////////////////////////////////////////////////////////
// GetLibrarySymbolAddress()
// returns a pointer to the data or function called "symbolname"
// from within the shared library called "libname"

void* GetLibrarySymbolAddress(const char* libname, const char* symbolname);


/////////////////////////////////////////////////////////////////
// GetLibraryHandle()
// opens, and returns the handle to, the library module
// called "libname"

LIBRARY_HANDLE GetLibraryHandle(const char* libname);


/////////////////////////////////////////////////////////////////
// GetHandleSymbolAddress()
// returns a pointer to the data or function called "symbolname"
// from within the shared library with handle "handle"

void* GetHandleSymbolAddress(LIBRARY_HANDLE handle, const char* symbolname);


/////////////////////////////////////////////////////////////////
// CloseLibrary()
//
// disassociates the given library handle from the calling process.

void CloseLibrary(LIBRARY_HANDLE);


////////////////////////////////////////////////////////////////
// SOError()
// returns the last error generated by one of the above functions.

const char* SOError();

LIBRARY_HANDLE findLib(std::string lib);
LIBRARY_HANDLE FindLibInPath(const std::string& lib, const std::string& path);

#endif
