#include <CCA/Components/Arches/Utility/UtilityFactory.h>
#include <CCA/Components/Arches/Utility/GridInfo.h>
#include <CCA/Components/Arches/Utility/TaskAlgebra.h>
#include <CCA/Components/Arches/Utility/SurfaceNormals.h>
#include <CCA/Components/Arches/Task/TaskInterface.h>

using namespace Uintah;

UtilityFactory::UtilityFactory()
{}

UtilityFactory::~UtilityFactory()
{}

void
UtilityFactory::register_all_tasks( ProblemSpecP& db )
{

  //GRID INFORMATION
  std::string tname = "grid_info";
  TaskInterface::TaskBuilder* tsk = scinew GridInfo::Builder( tname, 0 );
  register_task( tname, tsk );

  tname = "surface_normals";
  tsk = scinew SurfaceNormals::Builder( tname, 0 );
  register_task( tname, tsk );

  ProblemSpecP db_all_util = db->findBlock("Utilities");

  //<Utilities>
  if ( db_all_util ){
    for ( ProblemSpecP db_util = db_all_util->findBlock("utility"); db_util != 0;
          db_util = db_util->findNextBlock("utility")){

      std::string name;
      std::string type;
      db_util->getAttribute("label", name);
      db_util->getAttribute("type", type);

      if ( type == "variable_math" ){

        std::string grid_type="NA";
        if ( db_util->findBlock("grid").get_rep() != 0 ){
          db_util->findBlock("grid")->getAttribute("type", grid_type);
        } else {
          throw InvalidValue("Error: variable_math must specify a grid variable type", __FILE__, __LINE__);
        }

        TaskInterface::TaskBuilder* tsk;
        if ( grid_type == "CC" ){
          tsk = scinew TaskAlgebra<CCVariable<double> >::Builder(name, 0);
        } else if ( grid_type == "FX" ){
          tsk = scinew TaskAlgebra<SFCXVariable<double> >::Builder(name, 0);
        } else if ( grid_type == "FY" ){
          tsk = scinew TaskAlgebra<SFCYVariable<double> >::Builder(name, 0);
        } else if ( grid_type == "FZ" ){
          tsk = scinew TaskAlgebra<SFCZVariable<double> >::Builder(name, 0);
        } else {
          throw InvalidValue("Error: grid_type not recognized.",__FILE__,__LINE__);
        }
        register_task(name, tsk);

      } else {

        throw InvalidValue("Error: Utility type not recognized.",__FILE__,__LINE__);

      }

      assign_task_to_type_storage(name, type);

    }
  }
}

void
UtilityFactory::add_task( ProblemSpecP& db ){

  ProblemSpecP db_all_util = db->findBlock("Utilities");

  //<Utilities>
  if ( db_all_util ){
    for ( ProblemSpecP db_util = db_all_util->findBlock("utility"); db_util != 0;
          db_util = db_util->findNextBlock("utility")){

      std::string name;
      std::string type;
      db_util->getAttribute("label", name);
      db_util->getAttribute("type", type);

      if ( type == "variable_math" ){

        std::string grid_type="NA";
        if ( db_util->findBlock("grid").get_rep() != 0 ){
          db_util->findBlock("grid")->getAttribute("type", grid_type);
        } else {
          throw InvalidValue("Error: variable_math must specify a grid variable type", __FILE__, __LINE__);
        }

        TaskInterface::TaskBuilder* tsk;
        if ( grid_type == "CC" ){
          tsk = scinew TaskAlgebra<CCVariable<double> >::Builder(name, 0);
        } else if ( grid_type == "FX" ){
          tsk = scinew TaskAlgebra<SFCXVariable<double> >::Builder(name, 0);
        } else if ( grid_type == "FY" ){
          tsk = scinew TaskAlgebra<SFCYVariable<double> >::Builder(name, 0);
        } else if ( grid_type == "FZ" ){
          tsk = scinew TaskAlgebra<SFCZVariable<double> >::Builder(name, 0);
        } else {
          throw InvalidValue("Error: grid_type not recognized.",__FILE__,__LINE__);
        }
        register_task(name, tsk);

      } else {

        throw InvalidValue("Error: Utility type not recognized.",__FILE__,__LINE__);

      }

      assign_task_to_type_storage(name, type);

    }
  }
}

void
UtilityFactory::build_all_tasks( ProblemSpecP& db )
{

  // Grid spacing information
  TaskInterface* tsk = retrieve_task("grid_info");
  tsk->problemSetup(db);
  tsk->create_local_labels();

  tsk = retrieve_task("surface_normals");
  tsk->problemSetup(db);
  tsk->create_local_labels();

  //<Utilities>
  ProblemSpecP db_all_util = db->findBlock("Utilities");
  if ( db_all_util ){
    for ( ProblemSpecP db_util = db_all_util->findBlock("utility"); db_util != 0;
          db_util = db_util->findNextBlock("utility")){
      std::string name;
      db_util->getAttribute("label", name);
      TaskInterface* tsk = retrieve_task(name);
      tsk->problemSetup(db_util);
      tsk->create_local_labels();
    }
  }
}
