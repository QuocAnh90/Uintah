/*
 * The MIT License
 *
 * Copyright (c) 1997-2016 The University of Utah
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

/*
 *  uda2vis.cc: Provides an interface between VisIt and Uintah.
 *
 *  Written by:
 *   Department of Computer Science
 *   University of Utah
 *   April 2003-2007
 *
 */

#include <StandAlone/tools/uda2vis/uda2vis.h>

#include <Core/DataArchive/DataArchive.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <algorithm>

using namespace Uintah;


/////////////////////////////////////////////////////////////////////
// Utility functions for copying data from Uintah structures into
// simple arrays.
void copyIntVector(int to[3], const IntVector &from)
{
  to[0]=from[0];  to[1]=from[1];  to[2]=from[2];
}

void copyVector(double to[3], const Vector &from)
{
  to[0]=from[0];  to[1]=from[1];  to[2]=from[2];
}

void copyVector(double to[3], const Point &from)
{
  to[0]=from.x();  to[1]=from.y();  to[2]=from.z();
}

/////////////////////////////////////////////////////////////////////
// Utility functions for serializing Uintah data structures into
// a simple array for visit.
template <typename T>
int numComponents()
{
  return 1;
}

template <>
int numComponents<Vector>()
{
  return 3;
}

template <>
int numComponents<Stencil7>()
{
  return 7;
}

template <>
int numComponents<Stencil4>()
{
  return 4;
}

template <>
int numComponents<Point>()
{
  return 3;
}

template <>
int numComponents<Matrix3>()
{
  return 9;
}

template <typename T>
void copyComponents(double *dest, const T &src)
{
  (*dest) = (double)src;
}

template <>
void copyComponents<Vector>(double *dest, const Vector &src)
{
  dest[0] = (double)src[0];
  dest[1] = (double)src[1];
  dest[2] = (double)src[2];
}

template <>
void copyComponents<Stencil7>(double *dest, const Stencil7 &src)
{
  dest[0] = (double)src[0];
  dest[1] = (double)src[1];
  dest[2] = (double)src[2];
  dest[3] = (double)src[3];
  dest[4] = (double)src[4];
  dest[5] = (double)src[5];
  dest[6] = (double)src[6];
}

template <>
void copyComponents<Stencil4>(double *dest, const Stencil4 &src)
{
  dest[0] = (double)src[0];
  dest[1] = (double)src[1];
  dest[2] = (double)src[2];
  dest[3] = (double)src[3];
}

template <>
void copyComponents<Point>(double *dest, const Point &src)
{
  dest[0] = (double)src.x();
  dest[1] = (double)src.y();
  dest[2] = (double)src.z();
}

template <>
void copyComponents<Matrix3>(double *dest, const Matrix3 &src)
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      dest[i*3+j] = (double)src(i,j);
    }
  }
}


/////////////////////////////////////////////////////////////////////
// Open a data archive.
extern "C"
DataArchive* openDataArchive(const std::string& input_uda_name)
{
  DataArchive *archive = new DataArchive(input_uda_name);

  return archive;
}


/////////////////////////////////////////////////////////////////////
// Close a data archive - the visit plugin itself doesn't know about
// DataArchive::~DataArchive().
extern "C"
void closeDataArchive(DataArchive *archive)
{
  delete archive;
}


/////////////////////////////////////////////////////////////////////
// Get the grid for the current timestep, so we don't have to query
// it over and over.  We return a pointer to the GridP since the 
// visit plugin doesn't actually know about Grid's (or GridP's), and
// so the handle doesn't get destructed.
extern "C"
GridP* getGrid(DataArchive *archive, int timeStepNo)
{
  GridP *grid = new GridP(archive->queryGrid(timeStepNo));
  return grid;
}


/////////////////////////////////////////////////////////////////////
// Destruct the GridP, which will decrement the reference count.
extern "C"
void releaseGrid(GridP *grid)
{
  delete grid;
}


/////////////////////////////////////////////////////////////////////
// Get the time for each cycle.
extern "C"
std::vector<double> getCycleTimes(DataArchive *archive)
{

  // Get the times and indices.
  std::vector<int> index;
  std::vector<double> times;

  // query time info from dataarchive
  archive->queryTimesteps(index, times);

  return times;
} 


/////////////////////////////////////////////////////////////////////
// Get all the information that may be needed for the current timestep,
// including variable/material info, and level/patch info
extern "C"
TimeStepInfo* getTimeStepInfo(DataArchive *archive,
                              GridP *grid,
                              int timestep,
                              bool useExtraCells)
{
  int numLevels = (*grid)->numLevels();
  TimeStepInfo *stepInfo = new TimeStepInfo();
  stepInfo->levelInfo.resize(numLevels);

  // Get the variable information
  std::vector<std::string> vars;
  std::vector<const Uintah::TypeDescription*> types;
  archive->queryVariables(vars, types);
  stepInfo->varInfo.resize(vars.size());

  for (unsigned int i=0; i<vars.size(); i++) {
    VariableInfo &varInfo = stepInfo->varInfo[i];

    varInfo.name = vars[i];
    varInfo.type = types[i]->getName();

    // Query each level for material info until we find something
    for (int l=0; l<numLevels; l++) {
      LevelP level = (*grid)->getLevel(l);
      const Patch* patch = *(level->patchesBegin());
      ConsecutiveRangeSet matls =
        archive->queryMaterials(vars[i], patch, timestep);
      if (matls.size() > 0) {

        // Copy the list of materials
        for (ConsecutiveRangeSet::iterator matlIter = matls.begin();
             matlIter != matls.end(); matlIter++)
          varInfo.materials.push_back(*matlIter);

        // Don't query on any more levels
        break;
      }
    }
  }

  // Get level information
  for (int l=0; l<numLevels; l++) {
    LevelInfo &levelInfo = stepInfo->levelInfo[l];
    LevelP level = (*grid)->getLevel(l);

    copyIntVector(levelInfo.refinementRatio, level->getRefinementRatio());
    copyVector(levelInfo.spacing, level->dCell());
    copyVector(levelInfo.anchor, level->getAnchor());
    copyIntVector(levelInfo.periodic, level->getPeriodicBoundaries());

    // Patch info
    int numPatches = level->numPatches();
    levelInfo.patchInfo.resize(numPatches);

    for (int p=0; p<numPatches; p++) {
      const Patch* patch = level->getPatch(p);
      PatchInfo &patchInfo = levelInfo.patchInfo[p];

      // If the user wants to see extra cells, just include them and
      // let VisIt believe they are part of the original data. This is
      // accomplished by setting <meshtype>_low and <meshtype>_high to
      // the extra cell boundaries so that VisIt is none the wiser.
      if (useExtraCells)
      {
        patchInfo.setBounds(&patch->getExtraCellLowIndex()[0],
                            &patch->getExtraCellHighIndex()[0], "CC_Mesh");
        patchInfo.setBounds(&patch->getExtraNodeLowIndex()[0],
                            &patch->getExtraNodeHighIndex()[0], "NC_Mesh");
        patchInfo.setBounds(&patch->getExtraSFCXLowIndex()[0],
                            &patch->getExtraSFCXHighIndex()[0], "SFCX_Mesh");
        patchInfo.setBounds(&patch->getExtraSFCYLowIndex()[0],
                            &patch->getExtraSFCYHighIndex()[0], "SFCY_Mesh");
        patchInfo.setBounds(&patch->getExtraSFCZLowIndex()[0],
                            &patch->getExtraSFCZHighIndex()[0], "SFCZ_Mesh");
      }
      else
      {
        patchInfo.setBounds(&patch->getCellLowIndex()[0],
                            &patch->getCellHighIndex()[0], "CC_Mesh");
        patchInfo.setBounds(&patch->getNodeLowIndex()[0],
                            &patch->getNodeHighIndex()[0], "NC_Mesh");
        patchInfo.setBounds(&patch->getSFCXLowIndex()[0],
                            &patch->getSFCXHighIndex()[0], "SFCX_Mesh");
        patchInfo.setBounds(&patch->getSFCYLowIndex()[0],
                            &patch->getSFCYHighIndex()[0], "SFCY_Mesh");
        patchInfo.setBounds(&patch->getSFCZLowIndex()[0],
                            &patch->getSFCZHighIndex()[0], "SFCZ_Mesh");
      }

      // Get the processor id
      patchInfo.setProcId(archive->queryPatchwiseProcessor(patch, timestep));
    }
  }

  return stepInfo;
}


/////////////////////////////////////////////////////////////////////
// Read the grid data for the given index range
template<template <typename> class VAR, typename T>
static GridDataRaw* readGridData(DataArchive *archive,
                                 const Patch *patch,
                                 const LevelP level,
                                 std::string variable_name,
                                 int material,
                                 int timestep,
                                 int low[3],
                                 int high[3])
{
  IntVector ilow(low[0], low[1], low[2]);
  IntVector ihigh(high[0], high[1], high[2]);

  // this queries the entire patch, including extra cells and boundary cells
  VAR<T> var;
  archive->queryRegion(var, variable_name, material,
                       level.get_rep(), timestep, ilow, ihigh);

  //  IntVector low = var.getLowIndex();
  //  IntVector high = var.getHighIndex();

  GridDataRaw *gd = new GridDataRaw;
  gd->components = numComponents<T>();
  for (int i=0; i<3; i++) {
    gd->low[i] = low[i];
    gd->high[i] = high[i];
  }

  int n = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
  gd->data = new double[n*gd->components];

  T *p=var.getPointer();
  for (int i=0; i<n; i++)
    copyComponents<T>(&gd->data[i*gd->components], p[i]);
  
  return gd;
}


template<template<typename> class VAR>
GridDataRaw* getGridDataMainType(DataArchive *archive,
                                 const Patch *patch,
                                 const LevelP level,
                                 std::string variable_name,
                                 int material,
                                 int timestep,
                                 int low[3],
                                 int high[3],
                                 const Uintah::TypeDescription *subtype)
{
  switch (subtype->getType()) {
  case Uintah::TypeDescription::double_type:
    return readGridData<VAR, double>(archive, patch, level, variable_name,
                                     material, timestep, low, high);
  case Uintah::TypeDescription::float_type:
    return readGridData<VAR, float>(archive, patch, level, variable_name,
                                    material, timestep, low, high);
  case Uintah::TypeDescription::int_type:
    return readGridData<VAR, int>(archive, patch, level, variable_name,
                                  material, timestep, low, high);
  case Uintah::TypeDescription::Vector:
    return readGridData<VAR, Vector>(archive, patch, level, variable_name,
                                     material, timestep, low, high);
  case Uintah::TypeDescription::Stencil7:
    return readGridData<VAR, Stencil7>(archive, patch, level, variable_name,
                                       material, timestep, low, high);
  case Uintah::TypeDescription::Stencil4:
    return readGridData<VAR, Stencil4>(archive, patch, level, variable_name,
                                       material, timestep, low, high);
  case Uintah::TypeDescription::Matrix3:
    return readGridData<VAR, Matrix3>(archive, patch, level, variable_name,
                                      material, timestep, low, high);
  case Uintah::TypeDescription::bool_type:
  case Uintah::TypeDescription::short_int_type:
  case Uintah::TypeDescription::long_type:
  case Uintah::TypeDescription::long64_type:
    std::cerr << "Subtype " << subtype->getName() << " is not implemented..."
              << std::endl;
    return NULL;
  default:
    std::cerr << "Unknown subtype: "
          << subtype->getType() << "  "
          << subtype->getName() << std::endl;
    return NULL;
  }
}


extern "C"
GridDataRaw* getGridData(DataArchive *archive,
                         GridP *grid,
                         int level_i,
                         int patch_i,
                         std::string variable_name,
                         int material,
                         int timestep,
                         int low[3],
                         int high[3])
{
  LevelP level = (*grid)->getLevel(level_i);
  const Patch *patch = level->getPatch(patch_i);

  // figure out what the type of the variable we're querying is
  std::vector<std::string> vars;
  std::vector<const Uintah::TypeDescription*> types;
  archive->queryVariables(vars, types);

  const Uintah::TypeDescription* maintype = NULL;
  const Uintah::TypeDescription* subtype = NULL;

  for (unsigned int i=0; i<vars.size(); i++)
  {
    if (vars[i] == variable_name) {
      maintype = types[i];
      subtype = maintype->getSubType();
    }
  }

  if (!maintype || !subtype)
  {
    std::cerr << "couldn't find variable " << variable_name <<  std::endl;
    return NULL;
  }

  switch(maintype->getType())
  {
  case Uintah::TypeDescription::CCVariable:
    return getGridDataMainType<CCVariable>(archive, patch, level,
                                           variable_name, material, timestep,
                                           low, high, subtype);
  case Uintah::TypeDescription::NCVariable:
    return getGridDataMainType<NCVariable>(archive, patch, level,
                                           variable_name, material, timestep,
                                           low, high, subtype);
  case Uintah::TypeDescription::SFCXVariable:
    return getGridDataMainType<SFCXVariable>(archive, patch, level,
                                             variable_name, material, timestep
                                             , low, high, subtype);
  case Uintah::TypeDescription::SFCYVariable:
    return getGridDataMainType<SFCYVariable>(archive, patch, level,
                                             variable_name, material, timestep,
                                             low, high, subtype);
  case Uintah::TypeDescription::SFCZVariable:
    return getGridDataMainType<SFCZVariable>(archive, patch, level,
                                             variable_name, material, timestep,
                                             low, high, subtype);
  default:
    std::cerr << "Type is unknown." << std::endl;
    return NULL;
  }
}


/////////////////////////////////////////////////////////////////////
// Check to see if a variable exists.
extern "C"
bool variableExists(DataArchive *archive,
                    std::string variable_name)
{
  // figure out what the type of the variable we're querying is
  std::vector<std::string> vars;
  std::vector<const Uintah::TypeDescription*> types;
  archive->queryVariables(vars, types);

  const Uintah::TypeDescription* maintype = NULL;
  const Uintah::TypeDescription* subtype = NULL;

  for (unsigned int i=0; i<vars.size(); i++) {
    if (vars[i] == variable_name) {
      maintype = types[i];
      subtype = maintype->getSubType();
    }
  }

  return (maintype && subtype);
}


/////////////////////////////////////////////////////////////////////
// Read all the particle data for a given patch.
template<typename T>
ParticleDataRaw* readParticleData(DataArchive *archive,
                                  const Patch *patch,
                                  std::string variable_name,
                                  int material,
                                  int timestep)
{
  ParticleDataRaw *pd = new ParticleDataRaw;
  pd->components = numComponents<T>();
  pd->num = 0;

  // figure out which material we're interested in
  ConsecutiveRangeSet allMatls =
    archive->queryMaterials(variable_name, patch, timestep);

  ConsecutiveRangeSet matlsForVar;
  if (material<0) {
    matlsForVar = allMatls;
  }
  else {
    // make sure the patch has the variable - use empty material set
    // if it doesn't
    if (allMatls.size()>0 && allMatls.find(material) != allMatls.end())
      matlsForVar.addInOrder(material);
  }

  // first get all the particle subsets so that we know how many total
  // particles we'll have
  std::vector<ParticleVariable<T>*> particle_vars;
  for( ConsecutiveRangeSet::iterator matlIter =
         matlsForVar.begin(); matlIter != matlsForVar.end(); matlIter++ )
  {
    int matl = *matlIter;

    ParticleVariable<T> *var = new ParticleVariable<T>;
    archive->query(*var, variable_name, matl, patch, timestep);

    particle_vars.push_back(var);
    pd->num += var->getParticleSubset()->numParticles();
  }

  // copy all the data
  int pi = 0;
  pd->data = new double[pd->components * pd->num];
  for (unsigned int i=0; i<particle_vars.size(); i++)
  {
    ParticleSubset::iterator p;

    for (p = particle_vars[i]->getParticleSubset()->begin();
         p != particle_vars[i]->getParticleSubset()->end(); ++p)
    {
      //TODO: need to be able to read data as array of longs for
      //particle id, but copyComponents always reads double
      copyComponents<T>(&pd->data[pi*pd->components],
                        (*particle_vars[i])[*p]);
      ++pi;
    }
  }

  // cleanup
  for (unsigned int i=0; i<particle_vars.size(); i++)
    delete particle_vars[i];

  return pd;
}


extern "C"
ParticleDataRaw* getParticleData(DataArchive *archive,
                                 GridP *grid,
                                 int level_i,
                                 int patch_i,
                                 std::string variable_name,
                                 int material,
                                 int timestep)
{
  LevelP level = (*grid)->getLevel(level_i);
  const Patch *patch = level->getPatch(patch_i);

  // figure out what the type of the variable we're querying is
  std::vector<std::string> vars;
  std::vector<const Uintah::TypeDescription*> types;
  archive->queryVariables(vars, types);

  const Uintah::TypeDescription* maintype = NULL;
  const Uintah::TypeDescription* subtype = NULL;

  for (unsigned int i=0; i<vars.size(); i++) {
    if (vars[i] == variable_name) {
      maintype = types[i];
      subtype = maintype->getSubType();
    }
  }

  if (!maintype || !subtype) {
    std::cerr << "couldn't find variable " << variable_name <<  std::endl;
    return NULL;
  }

  switch (subtype->getType()) {
  case Uintah::TypeDescription::double_type:
    return readParticleData<double>(archive, patch, variable_name,
                                    material, timestep);
  case Uintah::TypeDescription::float_type:
    return readParticleData<float>(archive, patch, variable_name,
                                   material, timestep);
  case Uintah::TypeDescription::int_type:
    return readParticleData<int>(archive, patch, variable_name,
                                 material, timestep);
  case Uintah::TypeDescription::long64_type:
    return readParticleData<long64>(archive, patch, variable_name,
                                    material, timestep);
  case Uintah::TypeDescription::Point:
    return readParticleData<Point>(archive, patch, variable_name,
                                   material, timestep);
  case Uintah::TypeDescription::Vector:
    return readParticleData<Vector>(archive, patch, variable_name,
                                    material, timestep);
  case Uintah::TypeDescription::Stencil7:
    return readParticleData<Stencil7>(archive, patch, variable_name,
                                      material, timestep);
  case Uintah::TypeDescription::Stencil4:
    return readParticleData<Stencil4>(archive, patch, variable_name,
                                      material, timestep);
  case Uintah::TypeDescription::Matrix3:
    return readParticleData<Matrix3>(archive, patch, variable_name,
                                     material, timestep);
  default:
    std::cerr << "Unknown subtype for particle data: "
              << subtype->getName() << std::endl;
    return NULL;
  }
}

extern "C"
std::string getParticlePositionName(DataArchive *archive)
{
    return archive->getParticlePositionName();
}

/////////////////////////////////////////////////////////////////////
// Interface between VisIt's libsim and Uintah
namespace Uintah {

/////////////////////////////////////////////////////////////////////
// Get all the information that may be needed for the current timestep,
// including variable/material info, and level/patch info
TimeStepInfo* getTimeStepInfo2(SchedulerP schedulerP,
                               GridP gridP,
                               int timestep,
                               bool useExtraCells)
{
  DataWarehouse *dw = schedulerP->getLastDW();
  LoadBalancer  *lb = schedulerP->getLoadBalancer();

  int numLevels = gridP->numLevels();
  TimeStepInfo *stepInfo = new TimeStepInfo();
  stepInfo->levelInfo.resize(numLevels);

  // Get the variable information
  const std::set<const VarLabel*, VarLabel::Compare> varLabels =
    schedulerP->getInitialRequiredVars();

  // Get the material information
  Scheduler::VarLabelMaterialMap* pLabelMatlMap =
    schedulerP->makeVarLabelMaterialMap();

  stepInfo->varInfo.resize(varLabels.size());

  // Loop through all of the variables
  unsigned int i = 0;
  std::set<const VarLabel*, VarLabel::Compare>::iterator varIter;

  for (varIter = varLabels.begin(); varIter != varLabels.end(); ++varIter, ++i )
  {
    VariableInfo &varInfo = stepInfo->varInfo[i];
      
    varInfo.name = (*varIter)->getName();
    varInfo.type = (*varIter)->typeDescription()->getName();
    
    // Loop through all of the materials for this variable
    Scheduler::VarLabelMaterialMap::iterator matMapIter =
      pLabelMatlMap->find( varInfo.name );
    
    std::list< int > &materials = matMapIter->second;
    std::list< int >::iterator matIter;

    for (matIter = materials.begin(); matIter != materials.end(); ++matIter)
    {
      bool exists = false;

      // Check to make sure the variable exists on at least one patch
      // for at least one level.
      for (int l=0; l<numLevels; ++l)
      {
        LevelP level = gridP->getLevel(l);
        int numPatches = level->numPatches();

        for (int p=0; p<numPatches; ++p)
        {
          const Patch* patch = level->getPatch(p);

          if( dw->exists( *varIter, *matIter, patch ) )
          {
            // The variable exists on this level and patch.
            varInfo.materials.push_back( *matIter );
            exists = true;
            break;
          }
        }

        if( exists == true )
          break;
      }
    }
  }

  // Get the level information
  for (int l=0; l<numLevels; l++)
  {
    LevelInfo &levelInfo = stepInfo->levelInfo[l];
    LevelP level = gridP->getLevel(l);

    copyIntVector(levelInfo.refinementRatio, level->getRefinementRatio());
    copyVector(levelInfo.spacing, level->dCell());
    copyVector(levelInfo.anchor, level->getAnchor());
    copyIntVector(levelInfo.periodic, level->getPeriodicBoundaries());

    // Patch info
    int numPatches = level->numPatches();
    levelInfo.patchInfo.resize(numPatches);

    for (int p=0; p<numPatches; p++)
    {
      const Patch* patch = level->getPatch(p);
      PatchInfo &patchInfo = levelInfo.patchInfo[p];

      // If the user wants to see extra cells, just include them and
      // let VisIt believe they are part of the original data. This is
      // accomplished by setting <meshtype>_low and <meshtype>_high to
      // the extra cell boundaries so that VisIt is none the wiser.
      if (useExtraCells)
      {
        patchInfo.setBounds(&patch->getExtraCellLowIndex()[0],
                            &patch->getExtraCellHighIndex()[0], "CC_Mesh");
        patchInfo.setBounds(&patch->getExtraNodeLowIndex()[0],
                            &patch->getExtraNodeHighIndex()[0], "NC_Mesh");
        patchInfo.setBounds(&patch->getExtraSFCXLowIndex()[0],
                            &patch->getExtraSFCXHighIndex()[0], "SFCX_Mesh");
        patchInfo.setBounds(&patch->getExtraSFCYLowIndex()[0],
                            &patch->getExtraSFCYHighIndex()[0], "SFCY_Mesh");
        patchInfo.setBounds(&patch->getExtraSFCZLowIndex()[0],
                            &patch->getExtraSFCZHighIndex()[0], "SFCZ_Mesh");
      }
      else
      {
        patchInfo.setBounds(&patch->getCellLowIndex()[0],
                            &patch->getCellHighIndex()[0], "CC_Mesh");
        patchInfo.setBounds(&patch->getNodeLowIndex()[0],
                            &patch->getNodeHighIndex()[0], "NC_Mesh");
        patchInfo.setBounds(&patch->getSFCXLowIndex()[0],
                            &patch->getSFCXHighIndex()[0], "SFCX_Mesh");
        patchInfo.setBounds(&patch->getSFCYLowIndex()[0],
                            &patch->getSFCYHighIndex()[0], "SFCY_Mesh");
        patchInfo.setBounds(&patch->getSFCZLowIndex()[0],
                            &patch->getSFCZHighIndex()[0], "SFCZ_Mesh");
      }

      // Get the processor id
      patchInfo.setProcId( lb->getPatchwiseProcessorAssignment(patch) );
    }
  }

  return stepInfo;
}


// ****************************************************************************
//  Method: GetLevelAndLocalPatchNumber
//
//  Purpose:
//      Translates the global patch identifier to a refinement level and patch
//      number local to that refinement level.
//  
//  Programmer: sshankar, taken from implementation of the plugin, CHOMBO
//  Creation:   May 20, 2008
//
// ****************************************************************************
void GetLevelAndLocalPatchNumber(TimeStepInfo* stepInfo,
                                 int global_patch, 
                                 int &level, int &local_patch)
{
  int num_levels = stepInfo->levelInfo.size();
  int num_patches = 0;
  int tmp = global_patch;
  level = 0;

  while (level < num_levels)
  {
    num_patches = stepInfo->levelInfo[level].patchInfo.size();

    if (tmp < num_patches)
      break;

    tmp -= num_patches;
    level++;
  }

  local_patch = tmp;
}


// ****************************************************************************
//  Method: GetGlobalDomainNumber
//
//  Purpose:
//      Translates the level and local patch number into a global patch id.
//  
// ****************************************************************************
int GetGlobalDomainNumber(TimeStepInfo* stepInfo,
                          int level, int local_patch)
{
  int g = 0;

  for (int l=0; l<level; l++)
    g += stepInfo->levelInfo[l].patchInfo.size();
  g += local_patch;

  return g;
}


// ****************************************************************************
//  Method: getBounds
//
//  Purpose:
//   Returns the bounds for the given patch of the specified mesh 
//   based on periodicity and type.
//
//  Node centered data uses the same mesh as cell centered, 
//  but face centered meshes need an extra value for one axis,
//  unless they are periodic on that axis.
//
//  use patch_id=-1 to query all patches.
//
// ****************************************************************************
void getBounds(int low[3], int high[3],
               const std::string meshName,
               const LevelInfo &levelInfo,
               int patch_id)
{
  levelInfo.getBounds(low,high,meshName,patch_id);
  
  // debug5 << "getBounds(" << meshName << ",id=" << patch_id << ")=["
  //     << low[0] << "," << low[1] << "," << low[2] << "] to ["
  //     << high[0] << "," << high[1] << "," << high[2] << "]" << std::endl;
}


// ****************************************************************************
//  Method: CheckNaNs
//
//  Purpose:
//      Check for and warn about NaN values in the file.
//
//  Arguments:
//      num        data size
//      data       data
//      level      level that contains this patch
//      patch      patch that contains these cells
//
//  Returns:    none
//
//  Programmer: cchriste
//  Creation:   06.02.2012
//
//  Modifications:
const double NAN_REPLACE_VAL = 1.0E9;

void CheckNaNs(int num, double *data, int level, int patch)
{
  // replace nan's with a large negative number
  std::vector<int> nanCells;

  for (int i=0; i<num; i++) 
  {
    if (std::isnan(data[i]))
    {
      data[i] = NAN_REPLACE_VAL;
      nanCells.push_back(i);
    }
  }

  if (!nanCells.empty())
  {
    std::stringstream sstr;
    sstr << "NaNs exist in this file (patch " << patch
         << " of level " << level
         << "). They have been replaced by the value "
         <<  NAN_REPLACE_VAL << ".";

    if ((int)nanCells.size()>40)
    {
      sstr << std::endl << "First 20: ";

      for (int i=0;i<(int)nanCells.size() && i<20;i++)
        sstr << nanCells[i] << ",";

      sstr << std::endl << "Last 20: ";

      for (int i=(int)nanCells.size()-21;i<(int)nanCells.size();i++)
        sstr << nanCells[i] << ",";
    }
    else
    {
      for (int i=0;i<(int)nanCells.size();i++)
        sstr << nanCells[i] << ((int)nanCells.size()!=(i+1)?",":".");
    }

    std::cerr << "Uintah/VisIt Libsim warning : " << sstr.str() << std::endl;
    // avtCallback::IssueWarning(sstr.str().c_str());
  }
}


/////////////////////////////////////////////////////////////////////
// Read the grid data for the given index range
template<template <typename> class VAR, typename T>
static GridDataRaw* readGridData(SchedulerP schedulerP,
                                 const Patch *patch,
                                 const LevelP level,
                                 const VarLabel *varLabel,
                                 int material,
                                 int timestep,
                                 int low[3],
                                 int high[3])
{
  DataWarehouse *dw = schedulerP->getLastDW();

  std::string variable_name = varLabel->getName();

  IntVector ilow(low[0], low[1], low[2]);
  IntVector ihigh(high[0], high[1], high[2]);

  // this queries the entire patch, including extra cells and boundary cells
  VAR<T> var;

  schedulerP->getLoadBalancer()->getPatchwiseProcessorAssignment(patch);

  if( dw->exists( varLabel, material, patch ) )
    dw->getRegion( var, varLabel, material, level.get_rep(), ilow, ihigh );
  else
    return NULL;

  // IntVector low = var.getLowIndex();
  // IntVector high = var.getHighIndex();

  if( numComponents<T>() == 0 )
    return NULL;

  GridDataRaw *gd = new GridDataRaw;
  gd->components = numComponents<T>();

  for (int i=0; i<3; ++i)
  {
    gd->low[i] = low[i];
    gd->high[i] = high[i];
  }

  int n = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
  gd->data = new double[n*gd->components];

  const T *p=var.getPointer();
  for (int i=0; i<n; ++i)
    copyComponents<T>(&gd->data[i*gd->components], p[i]);
  
  return gd;
}


/////////////////////////////////////////////////////////////////////
// Read the grid data for a given patch.
template<template<typename> class VAR>
GridDataRaw* getGridDataMainType(SchedulerP schedulerP,
                                 const Patch *patch,
                                 const LevelP level,
                                 const VarLabel *varLabel,
                                 int material,
                                 int timestep,
                                 int low[3],
                                 int high[3],
                                 const Uintah::TypeDescription *subtype)
{
  switch (subtype->getType())
  {
  case Uintah::TypeDescription::double_type:
    return readGridData<VAR, double>(schedulerP, patch, level, varLabel,
                                     material, timestep, low, high);
  case Uintah::TypeDescription::float_type:
    return readGridData<VAR, float>(schedulerP, patch, level, varLabel,
                                    material, timestep, low, high);
  case Uintah::TypeDescription::int_type:
    return readGridData<VAR, int>(schedulerP, patch, level, varLabel,
                                  material, timestep, low, high);
  case Uintah::TypeDescription::Vector:
    return readGridData<VAR, Vector>(schedulerP, patch, level, varLabel,
                                     material, timestep, low, high);
  case Uintah::TypeDescription::Stencil7:
    return readGridData<VAR, Stencil7>(schedulerP, patch, level, varLabel,
                                       material, timestep, low, high);
  case Uintah::TypeDescription::Stencil4:
    return readGridData<VAR, Stencil4>(schedulerP, patch, level, varLabel,
                                       material, timestep, low, high);
  case Uintah::TypeDescription::Matrix3:
    return readGridData<VAR, Matrix3>(schedulerP, patch, level, varLabel,
                                      material, timestep, low, high);
  case Uintah::TypeDescription::bool_type:
  case Uintah::TypeDescription::short_int_type:
  case Uintah::TypeDescription::long_type:
  case Uintah::TypeDescription::long64_type:
    std::cerr << "Uintah/VisIt Libsim Error: "
              << "Subtype " << subtype->getName() << " is not implemented..."
              << std::endl;
    return NULL;
  default:
    std::cerr << "Uintah/VisIt Libsim Error: unknown subtype: "
              << subtype->getType() << "  for variable: "
              << subtype->getName() << std::endl;
    return NULL;
  }
}


/////////////////////////////////////////////////////////////////////
// Read the grid data for a given patch.
GridDataRaw* getGridData2(SchedulerP schedulerP,
                          GridP gridP,
                          int level_i,
                          int patch_i,
                          std::string variable_name,
                          int material,
                          int timestep,
                          int low[3],
                          int high[3])
{
  LevelP level = gridP->getLevel(level_i);
  const Patch *patch = level->getPatch(patch_i);

  // get the variable information
  const std::set<const VarLabel*, VarLabel::Compare> varLabels =
    schedulerP->getInitialRequiredVars();

  std::set<const VarLabel*, VarLabel::Compare>::iterator varIter;

  const VarLabel *varLabel;
  const Uintah::TypeDescription* maintype = NULL;
  const Uintah::TypeDescription* subtype = NULL;

  for (varIter = varLabels.begin(); varIter != varLabels.end(); varIter++)
  {
    if ((*varIter)->getName() == variable_name) {
    
      varLabel = (*varIter);    
      maintype = (*varIter)->typeDescription();
      subtype = (*varIter)->typeDescription()->getSubType();

      break;
    }
  }

  if (!maintype || !subtype) {
    std::cerr << "Uintah/VisIt Libsim Error: couldn't find variable "
              << variable_name <<  std::endl;
    return NULL;
  }

  switch(maintype->getType()) {
  case Uintah::TypeDescription::CCVariable:
    return getGridDataMainType<constCCVariable>(schedulerP, patch, level,
                                                varLabel, material, timestep,
                                                low, high, subtype);
  case Uintah::TypeDescription::NCVariable:
    return getGridDataMainType<constNCVariable>(schedulerP, patch, level,
                                                varLabel, material, timestep,
                                                low, high, subtype);
  case Uintah::TypeDescription::SFCXVariable:
    return getGridDataMainType<constSFCXVariable>(schedulerP, patch, level,
                                                  varLabel, material, timestep,
                                                  low, high, subtype);
  case Uintah::TypeDescription::SFCYVariable:
    return getGridDataMainType<constSFCYVariable>(schedulerP, patch, level,
                                                  varLabel, material, timestep,
                                                  low, high, subtype);
  case Uintah::TypeDescription::SFCZVariable:
    return getGridDataMainType<constSFCZVariable>(schedulerP, patch, level,
                                                  varLabel, material, timestep,
                                                  low, high, subtype);
  default:
    std::cerr << "Uintah/VisIt Libsim Error: unknown type: "
              << maintype->getName() << " for variable: "
              << variable_name << std::endl;
    return NULL;
  }
}


/////////////////////////////////////////////////////////////////////
// Read all the particle data for a given patch.
template<typename T>
ParticleDataRaw* readParticleData(SchedulerP schedulerP,
                                  const Patch *patch,
                                  const VarLabel *varLabel,
                                  int material,
                                  int timestep)
{
  DataWarehouse *dw = schedulerP->getLastDW();

  std::string variable_name = varLabel->getName();

  ParticleDataRaw *pd = new ParticleDataRaw;
  pd->components = numComponents<T>();
  pd->num = 0;

  // get the material information for all variables
  Scheduler::VarLabelMaterialMap* pLabelMatlMap =
    schedulerP->makeVarLabelMaterialMap();

  // get the materials for this variable
  Scheduler::VarLabelMaterialMap::iterator matMapIter =
    pLabelMatlMap->find( variable_name );

  // figure out which material we're interested in
  std::list< int > &allMatls = matMapIter->second;
  std::list< int > matlsForVar;

  if (material < 0)
  {
    matlsForVar = allMatls;
  }
  else
  {
    // make sure the patch has the variable - use empty material set
    // if it doesn't
    for (std::list< int >::iterator matIter = allMatls.begin();
         matIter != allMatls.end(); matIter++)
    {
      if( *matIter == material )
      {
        matlsForVar.push_back(material);
        break;
      }
    }
  }

  // first get all the particle subsets so that we know how many total
  // particles we'll have
  std::vector<constParticleVariable<T>*> particle_vars;

  for( std::list< int >::iterator matIter = matlsForVar.begin();
       matIter != matlsForVar.end(); matIter++ )
  {
    int matl = *matIter;

    constParticleVariable<T> *var = new constParticleVariable<T>;

    if( dw->exists( varLabel, matl, patch ) )
      dw->get( *var, varLabel, matl, patch);

    //archive->query(*var, variable_name, matl, patch, timestep);

    particle_vars.push_back(var);
    pd->num += var->getParticleSubset()->numParticles();
  }

  // figure out which material we're interested in
  // ConsecutiveRangeSet allMatls =
  //   archive->queryMaterials(variable_name, patch, timestep);

  // ConsecutiveRangeSet matlsForVar;

  // if (material < 0)
  // {
  //   matlsForVar = allMatls;
  // }
  // else
  // {
       // make sure the patch has the variable - use empty material set
       // if it doesn't
  //   if (0 < allMatls.size() && allMatls.find(material) != allMatls.end())
  //     matlsForVar.addInOrder(material);
  // }

  // first get all the particle subsets so that we know how many total
  // particles we'll have
  // std::vector<ParticleVariable<T>*> particle_vars;

  // for( ConsecutiveRangeSet::iterator matlIter = matlsForVar.begin();
  //      matlIter != matlsForVar.end(); matlIter++ )
  // {
  //   int matl = *matlIter;

  //   ParticleVariable<T> *var = new ParticleVariable<T>;
  //   archive->query(*var, variable_name, matl, patch, timestep);

  //   particle_vars.push_back(var);
  //   pd->num += var->getParticleSubset()->numParticles();
  // }

  // copy all the data
  int pi=0;
  pd->data = new double[pd->components * pd->num];

  for (unsigned int i=0; i<particle_vars.size(); i++)
  {
    ParticleSubset *pSubset = particle_vars[i]->getParticleSubset();

    for (ParticleSubset::iterator p = pSubset->begin();
         p != pSubset->end(); ++p)
    {

      //TODO: need to be able to read data as array of longs for
      //particle id, but copyComponents always reads double
      copyComponents<T>(&pd->data[pi*pd->components],
                        (*particle_vars[i])[*p]);
      pi++;
    }
  }

  // cleanup
  for (unsigned int i=0; i<particle_vars.size(); i++)
    delete particle_vars[i];

  return pd;
}


/////////////////////////////////////////////////////////////////////
// Read all the particle data for a given patch.
ParticleDataRaw* getParticleData2(SchedulerP schedulerP,
                                  GridP gridP,
                                  int level_i,
                                  int patch_i,
                                  std::string variable_name,
                                  int material,
                                  int timestep)
{
  LevelP level = gridP->getLevel(level_i);
  const Patch *patch = level->getPatch(patch_i);

  // get the variable information
  const std::set<const VarLabel*, VarLabel::Compare> varLabels =
    schedulerP->getInitialRequiredVars();

  std::set<const VarLabel*, VarLabel::Compare>::iterator varIter;

  const VarLabel *varLabel;
  const Uintah::TypeDescription* maintype = NULL;
  const Uintah::TypeDescription* subtype = NULL;
  

  for (varIter = varLabels.begin(); varIter != varLabels.end(); varIter++)
  {
    if ((*varIter)->getName() == variable_name) {

      varLabel = (*varIter);    
      maintype = (*varIter)->typeDescription();
      subtype = (*varIter)->typeDescription()->getSubType();

      break;
    }
  }

  if (!maintype || !subtype) {
    std::cerr << "Uintah/VisIt Libsim Error: couldn't find variable "
              << variable_name << std::endl;
    return NULL;
  }

  switch (subtype->getType()) {
  case Uintah::TypeDescription::double_type:
    return readParticleData<double>(schedulerP, patch, varLabel,
                                    material, timestep);
  case Uintah::TypeDescription::float_type:
    return readParticleData<float>(schedulerP, patch, varLabel,
                                   material, timestep);
  case Uintah::TypeDescription::int_type:
    return readParticleData<int>(schedulerP, patch, varLabel,
                                 material, timestep);
  case Uintah::TypeDescription::long64_type:
    return readParticleData<long64>(schedulerP, patch, varLabel,
                                    material, timestep);
  case Uintah::TypeDescription::Point:
    return readParticleData<Point>(schedulerP, patch, varLabel,
                                   material, timestep);
  case Uintah::TypeDescription::Vector:
    return readParticleData<Vector>(schedulerP, patch, varLabel,
                                    material, timestep);
  case Uintah::TypeDescription::Stencil7:
    return readParticleData<Stencil7>(schedulerP, patch, varLabel,
                                      material, timestep);
  case Uintah::TypeDescription::Stencil4:
    return readParticleData<Stencil4>(schedulerP, patch, varLabel,
                                      material, timestep);
  case Uintah::TypeDescription::Matrix3:
    return readParticleData<Matrix3>(schedulerP, patch, varLabel,
                                     material, timestep);
  default:
    std::cerr << "Uintah/VisIt Libsim Error: " 
              << "unknown subtype for particle data: " << subtype->getName()
              << " for vairable: " << variable_name << std::endl;
    return NULL;
  }
}

}
