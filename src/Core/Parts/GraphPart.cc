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


/*
 *  GraphPart.cc
 *
 *  Written by:
 *   Yarden Livnat
 *   Department of Computer Science
 *   University of Utah
 *   Sep 2001
 *
 *  Copyright (C) 2001 SCI Group
 */

#include <iostream>
#include <Core/Parts/GraphPart.h>

namespace SCIRun {
  
GraphPart::GraphPart( PartInterface *parent, const string &name)
  : Part( parent, name ), PartInterface( this, parent, "GraphGui" )
{
  interface_ = this;
  cerr <<"GP = " << this << endl;
}

GraphPart::~GraphPart()
{
}

void 
GraphPart::set_num_lines( int n )
{
  data_.resize(n);
}
  
void 
GraphPart::add_values( vector<double> &v)
{
  for (unsigned i=0; i<v.size(); i++)
    data_[i].push_back(v[i]);
  new_values( v );
};

} // namespace SCIRun


