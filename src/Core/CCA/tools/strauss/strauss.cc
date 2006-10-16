/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2004 Scientific Computing and Imaging Institute,
   University of Utah.

   
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

 

#include <algorithm>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>

#include <Core/Util/Environment.h>
#include <Core/Util/sci_system.h>
#include <Core/CCA/tools/strauss/strauss.h>

using namespace std;
using namespace SCIRun;

//**** CRCState
std::ostream& operator<<(CRCState& out, const Leader&) {
  out << out.leader;
  return out;
}
//**** End of CRCState

string absPath(string& file)
{
  if(file.substr(0,1) != "/") {
    char buf[100];
    string absfile = std::string(getcwd(buf,100))+"/"+file;
    return absfile;
  }
  else return file;
}

Strauss::Strauss(string plugin, string hdrplugin, string portspec, 
		 string header, string implementation,
		 string util, string templateArgv)
  :header(header), implementation(implementation), plugin(plugin), 
   hdrplugin(hdrplugin), portspec(portspec), util(util), templateArgv(templateArgv)
{
  emitted = false;

  //Make paths absolute for plugin and portSpec 
  plugin = absPath(plugin);
  portspec = absPath(portspec);

  //Create name the bridge component class
  ostringstream o;
  srand(time(NULL));
  o << "Bridge" << 1+(int)(100000.0*rand()/(RAND_MAX+1.0));
  bridgeComponent = o.str();
}

Strauss::~Strauss()
{
}

int Strauss::emit()
{
  string scratchfile = "scratchfile";
  int bufsize=512;
  char *buf=new char[bufsize];

  if(emitted) return -1;
 
  impl << "//This file was automatically generated by Strauss. Do not edit directly!!\n";
  impl << "\n";
  impl << "#include \"" << header << "\"\n";
 
  hdr << "//This file was automatically generated by Strauss. Do not edit directly!!\n";
  hdr << "\n";
  hdr << "#ifndef STRAUSS_GENERATED_" << bridgeComponent << "\n";
  hdr << "#define STRAUSS_GENERATED_" << bridgeComponent << "\n";

  //create the name of the make component function
  string makename = header.substr(0,header.rfind("."));
  makename = makename.substr(header.rfind("/")+1);

  //Run SUPA
  int status;
  string execline;

  string srcdir = sci_getenv("SCIRUN_SRCDIR");
  string executable(srcdir + string("/Core/CCA/tools/scim/scim"));

  ///////////////////////////
  //IMPLEMENTATION:
  execline = executable + " " + "-t " + plugin + " -o " + scratchfile + " -r " + bridgeComponent + " -T " + makename;
  if(util != "") execline += " -R " + util;
  execline += " " + portspec;
  cerr << execline << "\n";
  status = sci_system(execline.c_str());
  if(status!=0) {
    return -1;
  }
  
  ifstream ifile(scratchfile.c_str());
  while(!ifile.eof()) {
    ifile.getline(buf,bufsize-1);
    string sbuf(buf);
    impl << sbuf << "\n";
  }
 

  ///////////////////////
  //HEADER
  execline = executable + " " + "-t " + hdrplugin + " -o " + scratchfile + " -r " + bridgeComponent; 
  if(util != "") execline += " -R " + util;
  execline += " " + portspec;
  cerr << execline << "\n";
  status = sci_system(execline.c_str());
  if(status!=0) {
    return -1;
  }
  
  ifstream hfile(scratchfile.c_str());
  while(!hfile.eof()) {
    hfile.getline(buf,bufsize-1);
    string sbuf(buf);
    hdr << sbuf << "\n";
  }
  ////////////////////

  delete buf;
  hdr << "#endif\n";

  /////
  // Calculate CRC and commit to files:
  impl.calcCRC(bridgeComponent); 
  hdr.calcCRC(bridgeComponent);
  emitted = true;
  commitToFiles();

  return 0;
}

unsigned long Strauss::getImplCRC()
{
  return impl.crc;
}

unsigned long Strauss::getHdrCRC()
{
  return hdr.crc;
}

//PRIVATE:

void Strauss::commitToFiles()
{
  if(emitted) {
    ofstream fHeader;
    ofstream fImpl;
    fHeader.open(header.c_str()); 
    fImpl.open(implementation.c_str());
    if((!fHeader)||(!fImpl)) {
      cerr << "ERROR in: ./strauss/main.cc: Couldn't open output file\n";
      exit(1);
    }

    fImpl << "//*****CRC=" << impl.crc << "\n";
    fImpl << impl.str();
    fHeader << "//*****CRC=" << hdr.crc << "\n";
    fHeader << hdr.str();
  }
}


