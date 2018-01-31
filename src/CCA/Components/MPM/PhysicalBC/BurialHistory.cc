/*
 * The MIT License
 *
 * Copyright (c) 1997-2018 The University of Utah
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

#include <vector>
#include <CCA/Components/MPM/PhysicalBC/BurialHistory.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Uintah;
using namespace std;

BurialHistory::BurialHistory(/*const ProcessorGroup* myworld*/)
{
}

BurialHistory::~BurialHistory()
{ 
}

void BurialHistory::populate(ProblemSpecP& ps)
{
  ProblemSpecP root = ps->getRootNode();
  ProblemSpecP burHist = root->findBlock("BurialHistory");
/*
  if (!burHist) 
    throw ProblemSetupException("**ERROR** No burial history specified.", 
                                 __FILE__, __LINE__);
*/
  if (burHist){
   burHist->getWithDefault("pressure_conversion_factor", 
                          d_pressure_conversion_factor, 1.0);
   for( ProblemSpecP timePoint = burHist->findBlock("time_point");
       timePoint != nullptr;
       timePoint = timePoint->findNextBlock("time_point") ) {
     double time      = 0.0;
     double depth     = 9.9e99;
     double temp      = 0.0;
     double fluidP    = 0.0;
     double fluidOP   = 0.0;
     double effStress = 0.0;
     double waterSat  = 0.0;
     timePoint->require("time_Ma",               time);
     timePoint->require("depth_m",               depth);
     timePoint->require("temperature_K",         temp);
     timePoint->require("fluidOverPressure_bar", fluidOP);
     timePoint->require("fluidPressure_bar",     fluidP);
     timePoint->require("effectiveStress_bar",   effStress);
     timePoint->require("waterSaturation_pct",   waterSat);

     d_time_Ma.push_back(time);
     d_depth_m.push_back(depth);
     d_temperature_K.push_back(temp);
     d_fluidOverPressure_bar.push_back(fluidOP);
     d_fluidPressure_bar.push_back(fluidP);
     d_effectiveStress_bar.push_back(effStress);
     d_waterSaturation_pct.push_back(waterSat);
   }
   burHist->getWithDefault("current_index", d_CI, d_time_Ma.size()-1);
  } // endif
}

void BurialHistory::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP lc_ps = ps->appendChild("BurialHistory");
  lc_ps->appendElement("pressure_conversion_factor",
                      d_pressure_conversion_factor);
  lc_ps->appendElement("current_index", d_CI);
  for (int i = 0; i<(int)d_time_Ma.size();i++) {
    ProblemSpecP time_ps = lc_ps->appendChild("time_point");
    time_ps->appendElement("time_Ma",              d_time_Ma[i]);
    time_ps->appendElement("depth_m",              d_depth_m[i]);
    time_ps->appendElement("temperature_K",        d_temperature_K[i]);
    time_ps->appendElement("fluidPressure_bar",    d_fluidPressure_bar[i]);
    time_ps->appendElement("fluidOverPressure_bar",d_fluidOverPressure_bar[i]);
    time_ps->appendElement("effectiveStress_bar",  d_effectiveStress_bar[i]);
    time_ps->appendElement("waterSaturation_pct",  d_waterSaturation_pct[i]);
  }
}
