/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2004 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/



/*
 *  AssertionFailed.h: Generic exception for internal errors
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   July 1999
 *
 *  Copyright (C) 1999 SCI Group
 */

#include <Core/Exceptions/AssertionFailed.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

namespace SCIRun {

AssertionFailed::AssertionFailed(const char* message,
				 const char* file,
				 int lineno)
{
    size_t len = strlen(message)+strlen(file)+strlen(message_)+100;
    message_ = (char*)malloc(len);
    sprintf(message_, "%s (file: %s, line: %d)\n%s", message, file, lineno, stacktrace_);
}

AssertionFailed::AssertionFailed(const AssertionFailed& copy)
    : message_(strdup(copy.message_))
{
}

AssertionFailed::~AssertionFailed()
{
    free(message_);
}

const char* AssertionFailed::message() const
{
    return message_;
}

const char* AssertionFailed::type() const
{
    return "AssertionFailed";
}

} // End namespace SCIRun
