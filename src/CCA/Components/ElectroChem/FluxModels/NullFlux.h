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

#ifndef CCA_COMPONENTS_ELECTROCHEM_FLUXMODELS_NULLFLUX_H
#define CCA_COMPONENTS_ELECTROCHEM_FLUXMODELS_NULLFLUX_H

#include <CCA/Components/ElectroChem/FluxModels/FluxModel.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

using namespace Uintah;

namespace ElectroChem {
  class NullFlux : public FluxModel {
    MaterialSet*    d_one_mat_set;
    MaterialSubset* d_one_mat_subset;
    public:
      NullFlux(ProblemSpecP& ps);
      virtual ~NullFlux();

      virtual void AddComputesAndRequires(const PatchSet* patches,
                                                Task*     task);
      virtual void ComputeFlux(const PatchSubset*   patches,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw);
  }; // End class NullFlux
} // End namspace ElectroChem
#endif // End CCA_COMPONENTS_ELECTROCHEM_FLUXMODELS_NULLFLUX_H
