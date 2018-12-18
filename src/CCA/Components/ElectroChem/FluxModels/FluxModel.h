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

#ifndef CCA_COMPONENTS_ELECTROCHEM_FLUXMODELS_FLUXMODEL_H
#define CCA_COMPONENTS_ELECTROCHEM_FLUXMODELS_FLUXMODEL_H

#include <CCA/Components/ElectroChem/ECLabel.h>

#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ComputeSet.h>

using namespace Uintah;

namespace ElectroChem {
  class FluxModel {
    MaterialSet*    d_one_mat_set    {nullptr};
    MaterialSubset* d_one_mat_subset {nullptr};
    ECLabel*        d_eclabel        {nullptr};

    public:
      FluxModel() { }
      virtual ~FluxModel() { }

      void SetMaterialSets(MaterialSet* set, MaterialSubset* subset){
        d_one_mat_set    = set;
        d_one_mat_subset = subset;
      }

      void SetECLabel(ECLabel* ec_label) { d_eclabel = ec_label; }

      virtual void AddComputesAndRequires(const PatchSet* patches,
                                                Task*     task) = 0;

      virtual void ComputeFlux(const PatchSubset*   patches,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw) = 0;
  }; // End class BasicFlux
} // End namspace ElectroChem
#endif // End CCA_COMPONENTS_ELECTROCHEM_FLUXMODELS_FLUXMODEL_H
