#ifndef Wasatch_TimeStepper_h
#define Wasatch_TimeStepper_h

#include <spatialops/structured/FVStaggeredTypes.h>

#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/ComputeSet.h>

#include <expression/Expr_ExpressionID.h>
#include <expression/FieldManager.h> // field type conversion tools
#include <expression/ExpressionFactory.h>
#include <expression/PlaceHolderExpr.h>

#include "PatchInfo.h"
#include "FieldAdaptor.h"

namespace Uintah{
  class ProcessorGroup;
  class DataWarehouse;
}

namespace Wasatch{

  /**
   *  \class  TimeStepper
   *  \author James C. Sutherland
   *  \date   June 2010
   *
   *  \brief Support for integrating a set of transport equations
   *         (explicit time integration methods for now).
   */
  class TimeStepper
  {

    /**
     *  \struct FieldInfo
     *  \author James C. Sutherland
     *
     *  \brief provides strongly typed information about a field.
     *         These are used to provide information about what fields
     *         we are advancing in the time integrator, and their
     *         associated RHS expressions.
     */
    template<typename FieldT>
    struct FieldInfo
    {
      std::string varname;
      Uintah::VarLabel* varLabel;
      Uintah::VarLabel* rhsLabel;
      FieldInfo( const std::string& name,
                 Uintah::VarLabel* const vl,
                 Uintah::VarLabel* const rhsl )
        : varname( name ), varLabel( vl ), rhsLabel( rhsl )
      {}
    };

    typedef std::vector< FieldInfo<SpatialOps::structured::SVolField> > ScalarFields;
    typedef std::vector< FieldInfo<SpatialOps::structured::XVolField> > XVolFields;
    typedef std::vector< FieldInfo<SpatialOps::structured::YVolField> > YVolFields;
    typedef std::vector< FieldInfo<SpatialOps::structured::ZVolField> > ZVolFields;

    ScalarFields scalarFields_;  ///< A vector of the scalar fields being solved by this time integrator.
    XVolFields   xVolFields_;    ///< A vector of the x-volume fields being solved by this time integrator.
    YVolFields   yVolFields_;    ///< A vector of the y-volume fields being solved by this time integrator.
    ZVolFields   zVolFields_;    ///< A vector of the z-volume fields being solved by this time integrator.

    typedef std::vector< Expr::ExpressionID > RHSIDList;
    RHSIDList rhsIDs_;  ///< A list of all of the RHS evaluators associated with this integrator.

    Expr::ExpressionFactory* const factory_;  ///< the factory that is associated with this time stepper.
    const Uintah::VarLabel* const deltaTLabel_;  ///< label for the time step variable.

    /**
     *  \brief used internally to obtain the appropriate vector
     *         (e.g. scalarFields_) given the type of field we are
     *         considering.
     */
    template<typename FieldT>
    std::vector< FieldInfo<FieldT> >& field_info_selctor();

    /**
     *  \brief the call-back for Uintah to execute this.
     */
    void update_variables( const Uintah::ProcessorGroup* const,
                           const Uintah::PatchSubset* const,
                           const Uintah::MaterialSubset* const,
                           Uintah::DataWarehouse* const,
                           Uintah::DataWarehouse* const );

  public:

    /**
     *  \brief Construct a TimeStepper object to advance equations forward in time
     *
     *  \param factory the ExpressionFactory that will be used to
     *         construc the trees for any transport equations added to
     *         this library.  The same factory should be used when
     *         constructing the expressions in each transport
     *         equation.
     */
    TimeStepper( const Uintah::VarLabel* deltaTLabel,
                 Expr::ExpressionFactory& factory );

    /**
     *  \brief Add a transport equation to this TimeStepper
     *
     *  \param solnVarName the name of the solution variable for this transport equation.
     *
     *  \param rhsID the Expr::ExpressionID for the right-hand-side of this transport equation.
     *
     *  This method is strongly typed to ensure that the solution
     *  variables are advanced properly and to guarantee compatibility
     *  with the Expression library.
     *
     *  The TimeStepper maintains an ExpressionTree object that
     *  calculates the RHS for all transport equations that are
     *  managed by this TimeStepper.  Adding an equation augments this
     *  tree and causes a recompilation of it.
     */
    template<typename FieldT>
    inline void add_equation( const std::string& solnVarName,
                              Expr::ExpressionID rhsID );

    /**
     *  \brief schedule the tasks associated with this TimeStepper
     */
    void create_tasks( const Expr::ExpressionID timeID,
                       const PatchInfoMap&,
                       const Uintah::PatchSet* const localPatches,
                       const Uintah::MaterialSet* const materials,
                       Uintah::SchedulerP& sched );
  };

  //------------------------------------------------------------------

  template<>
  inline std::vector< TimeStepper::FieldInfo<SpatialOps::structured::SVolField> >&
  TimeStepper::field_info_selctor<SpatialOps::structured::SVolField>()
  {
    return scalarFields_;
  }
  template<>
  inline std::vector<TimeStepper::FieldInfo<SpatialOps::structured::XVolField> >&
  TimeStepper::field_info_selctor<SpatialOps::structured::XVolField>()
  {
    return xVolFields_;
  }
  template<>
  inline std::vector<TimeStepper::FieldInfo<SpatialOps::structured::YVolField> >&
  TimeStepper::field_info_selctor<SpatialOps::structured::YVolField>()
  {
    return yVolFields_;
  }
  template<>
  inline std::vector<TimeStepper::FieldInfo<SpatialOps::structured::ZVolField> >&
  TimeStepper::field_info_selctor<SpatialOps::structured::ZVolField>()
  {
    return zVolFields_;
  }

  //------------------------------------------------------------------

  template<typename FieldT>
  void
  TimeStepper::add_equation( const std::string& solnVarName,
                             Expr::ExpressionID rhsID )
  {
    const std::string& rhsName = factory_->get_registry().get_label(rhsID).name();
    const Uintah::TypeDescription* typeDesc = getUintahFieldTypeDescriptor<FieldT>();
    const Uintah::IntVector ghostDesc       = getUintahGhostDescriptor<FieldT>();
    Uintah::VarLabel* solnVarLabel = Uintah::VarLabel::create( solnVarName, typeDesc, ghostDesc );
    Uintah::VarLabel* rhsVarLabel  = Uintah::VarLabel::create( rhsName,     typeDesc, ghostDesc );
    std::vector< FieldInfo<FieldT> >& fields = field_info_selctor<FieldT>();
    fields.push_back( FieldInfo<FieldT>( solnVarName, solnVarLabel, rhsVarLabel ) );
    rhsIDs_.push_back( rhsID );

    typedef Expr::PlaceHolder<FieldT>  FieldExpr;
    factory_->register_expression( Expr::Tag(solnVarName,Expr::STATE_N  ), new typename FieldExpr::Builder() );
    factory_->register_expression( Expr::Tag(solnVarName,Expr::STATE_NP1), new typename FieldExpr::Builder() );
  }

  //==================================================================

} // namespace Wasatch

#endif // Wasatch_TimeStepper_h
