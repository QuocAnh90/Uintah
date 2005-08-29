
#include <Packages/Uintah/CCA/Components/ICE/Thermo/ConstantThermo.h>
#include <Packages/Uintah/Core/Grid/Variables/CellIterator.h>

using namespace Uintah;

ConstantThermo::ConstantThermo(ProblemSpecP& ps, ModelSetup*, ICEMaterial* ice_matl)
  : ThermoInterface(ice_matl)
{
  ps->require("thermal_conductivity",d_thermalConductivity);
  ps->require("specific_heat",d_specificHeat);
  ps->require("gamma",d_gamma);
}

ConstantThermo::~ConstantThermo()
{
}

void ConstantThermo::scheduleInitializeThermo(SchedulerP& sched,
                                              const PatchSet* patches)
{
  // No initialization
}

void ConstantThermo::scheduleReactions(SchedulerP& sched,
                                       const PatchSet* patches)
{
  // No reactions
}

void ConstantThermo::addTaskDependencies_thermalDiffusivity(Task* t, State state,
                                                            int numGhostCells)
{
  // No additional requirements
}

void ConstantThermo::addTaskDependencies_thermalConductivity(Task* t, State state,
                                                             int numGhostCells)
{
  // No additional requirements
}

void ConstantThermo::addTaskDependencies_cp(Task* t, State state,
                                            int numGhostCells)
{
  // No additional requirements
}

void ConstantThermo::addTaskDependencies_cv(Task* t, State state,
                                            int numGhostCells)
{
  // No additional requirements
}

void ConstantThermo::addTaskDependencies_gamma(Task* t, State state,
                                               int numGhostCells)
{
  // No additional requirements
}

void ConstantThermo::addTaskDependencies_R(Task* t, State state,
                                           int numGhostCells)
{
  // No additional requirements
}

void ConstantThermo::addTaskDependencies_Temp(Task* t, State state,
                                              int numGhostCells)
{
  // No additional requirements
}

void ConstantThermo::addTaskDependencies_int_eng(Task* t, State state,
                                                 int numGhostCells)
{
  // No additional requirements
}

void ConstantThermo::compute_thermalDiffusivity(CellIterator iter,
                                                CCVariable<double>& thermalDiffusivity,
                                                DataWarehouse* dw, const Patch* patch,
                                                int matl, int numGhostCells,
                                                constCCVariable<double>& int_eng,
                                                constCCVariable<double>& sp_vol)
{
  double cp = d_specificHeat * d_gamma;
  double factor = d_thermalConductivity/cp;
  for(;!iter.done();iter++)
    thermalDiffusivity[*iter] = factor * sp_vol[*iter];
}

void ConstantThermo::compute_thermalConductivity(CellIterator iter,
                                                 CCVariable<double>& thermalConductivity,
                                                 DataWarehouse* dw, const Patch* patch,
                                                 int matl, int numGhostCells,
                                                 constCCVariable<double>& int_eng)
{
  for(;!iter.done();iter++)
    thermalConductivity[*iter] = d_thermalConductivity;
}

void ConstantThermo::compute_cp(CellIterator iter, CCVariable<double>& cp,
                                DataWarehouse* dw, const Patch* patch,
                                int matl, int numGhostCells,
                                constCCVariable<double>& int_eng)
{
  double tmp = d_specificHeat * d_gamma;
  for(;!iter.done();iter++)
    cp[*iter] = tmp;
}

void ConstantThermo::compute_cv(CellIterator iter, CCVariable<double>& cv,
                                DataWarehouse* dw, const Patch* patch,
                                int matl, int numGhostCells,
                                constCCVariable<double>&)
{
  for(;!iter.done();iter++)
    cv[*iter] = d_specificHeat;
}

void ConstantThermo::compute_gamma(CellIterator iter, CCVariable<double>& gamma,
                                   DataWarehouse* dw, const Patch* patch,
                                   int matl, int numGhostCells,
                                   constCCVariable<double>&)
{
  for(;!iter.done();iter++)
    gamma[*iter] = d_gamma;
}

void ConstantThermo::compute_R(CellIterator iter, CCVariable<double>& R,
                               DataWarehouse* dw, const Patch* patch,
                               int matl, int numGhostCells,
                               constCCVariable<double>&)
{
  double tmp = (d_gamma-1) * d_specificHeat;
  for(;!iter.done();iter++)
    R[*iter] = tmp;
}

void ConstantThermo::compute_Temp(CellIterator iter, CCVariable<double>& temp,
                                  DataWarehouse* dw, const Patch* patch,
                                  int matl, int numGhostCells,
                                  constCCVariable<double>& int_eng)
{
  double factor = 1./d_specificHeat;
  for(;!iter.done();iter++)
    temp[*iter] = int_eng[*iter] * factor;
}

void ConstantThermo::compute_int_eng(CellIterator iter, CCVariable<double>& int_eng,
                                     DataWarehouse* dw, const Patch* patch,
                                     int matl, int numGhostCells,
                                     constCCVariable<double>& temp,
                                     constCCVariable<double>&)
{
  for(;!iter.done();iter++)
    int_eng[*iter] = temp[*iter] * d_specificHeat;
}
