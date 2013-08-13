/*
 * The MIT License
 *
 * Copyright (c) 2012 The University of Utah
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
#include "MomentumTransportEquation.h"

// -- Uintah includes --//
#include <CCA/Ports/SolverInterface.h>
#include <Core/Exceptions/ProblemSetupException.h>

//-- Wasatch includes --//
#include <CCA/Components/Wasatch/TagNames.h>
#include <CCA/Components/Wasatch/Operators/OperatorTypes.h>
#include <CCA/Components/Wasatch/Expressions/TimeDerivative.h>
#include <CCA/Components/Wasatch/Expressions/MomentumPartialRHS.h>
#include <CCA/Components/Wasatch/Expressions/MomentumRHS.h>
#include <CCA/Components/Wasatch/Expressions/Strain.h>
#include <CCA/Components/Wasatch/Expressions/Dilatation.h>
#include <CCA/Components/Wasatch/Expressions/Turbulence/TurbulentViscosity.h>
#include <CCA/Components/Wasatch/Expressions/Turbulence/StrainTensorBase.h>
#include <CCA/Components/Wasatch/Expressions/Turbulence/StrainTensorMagnitude.h>
#include <CCA/Components/Wasatch/Expressions/Turbulence/DynamicSmagorinskyCoefficient.h>
#include <CCA/Components/Wasatch/Expressions/MMS/Functions.h>
#include <CCA/Components/Wasatch/Expressions/BoundaryConditions/BoundaryConditionBase.h>
#include <CCA/Components/Wasatch/Expressions/BoundaryConditions/BCCopier.h>

#include <CCA/Components/Wasatch/Expressions/EmbeddedGeometry/EmbeddedGeometryHelper.h>
#include <CCA/Components/Wasatch/Expressions/PrimVar.h>
#include <CCA/Components/Wasatch/Expressions/PressureSource.h>
#include <CCA/Components/Wasatch/Expressions/VelEst.h>
#include <CCA/Components/Wasatch/Expressions/WeakConvectiveTerm.h>
#include <CCA/Components/Wasatch/Expressions/ExprAlgebra.h>
#include <CCA/Components/Wasatch/Expressions/PostProcessing/InterpolateExpression.h>
#include <CCA/Components/Wasatch/Expressions/ConvectiveFlux.h>
#include <CCA/Components/Wasatch/Expressions/Pressure.h>
#include <CCA/Components/Wasatch/ConvectiveInterpolationMethods.h>
#include <CCA/Components/Wasatch/OldVariable.h>
#include <CCA/Components/Wasatch/FieldTypes.h>
#include <CCA/Components/Wasatch/ParseTools.h>

#include <CCA/Components/Wasatch/ReductionHelper.h>

#include <CCA/Components/Wasatch/Expressions/PostProcessing/KineticEnergy.h>

//-- ExprLib Includes --//
#include <expression/ExprLib.h>

using std::string;

namespace Wasatch{

  //==================================================================
  
  void register_turbulence_expressions (const TurbulenceParameters& turbParams,
                                        Expr::ExpressionFactory& factory,
                                        const Expr::TagList& velTags,
                                        const Expr::Tag densTag,
                                        const bool isConstDensity) {

    const TagNames& tagNames = TagNames::self();
    
    Expr::Tag strTsrMagTag      = Expr::Tag();
    Expr::Tag waleTsrMagTag     = Expr::Tag();
    Expr::Tag dynSmagCoefTag    = Expr::Tag();
    Expr::Tag vremanTsrMagTag   = Expr::Tag();
    const Expr::Tag turbViscTag = tagNames.turbulentviscosity;

    // Disallow users from using turbulence models in 1 or 2 dimensions
    if (!( velTags[0]!=Expr::Tag() && velTags[1]!=Expr::Tag() && velTags[2]!=Expr::Tag() )) {
      std::ostringstream msg;
      msg << "ERROR: You cannot use a turbulence model in one or two dimensions. Please revise your input file and make sure that you specify all three velocity/momentum components." << std::endl;
      throw Uintah::ProblemSetupException( msg.str(), __FILE__, __LINE__ );
    }

    // we have turbulence turned on. create an expression for the strain tensor magnitude. this is used by all eddy viscosity models
    switch (turbParams.turbModelName) {
        
        // ---------------------------------------------------------------------
      case SMAGORINSKY: {
        strTsrMagTag = tagNames.straintensormag;//( "StrainTensorMagnitude", Expr::STATE_NONE );
        if( !factory.have_entry( strTsrMagTag ) ){
          typedef StrainTensorSquare::Builder StrTsrMagT;
          factory.register_expression( scinew StrTsrMagT(strTsrMagTag,
                                                         tagNames.tauxx,tagNames.tauyx,tagNames.tauzx,
                                                         tagNames.tauyy,tagNames.tauzy,
                                                         tagNames.tauzz) );
        }
      }
        break;

        // ---------------------------------------------------------------------
      case VREMAN: {
        vremanTsrMagTag = tagNames.vremantensormag;
        if( !factory.have_entry( vremanTsrMagTag ) ){
          typedef VremanTensorMagnitude::Builder VremanTsrMagT;
          factory.register_expression( scinew VremanTsrMagT(vremanTsrMagTag, velTags ) );
        }
      }
        break;
        
        // ---------------------------------------------------------------------
      case WALE: {
        strTsrMagTag = tagNames.straintensormag;
        if( !factory.have_entry( strTsrMagTag ) ){
          typedef StrainTensorSquare::Builder StrTsrMagT;
          factory.register_expression( scinew StrTsrMagT(strTsrMagTag,
                                                         tagNames.tauxx,tagNames.tauyx,tagNames.tauzx,
                                                         tagNames.tauyy,tagNames.tauzy,
                                                         tagNames.tauzz) );
        }
        
        // if WALE model is turned on, then create an expression for the square velocity gradient tensor
        waleTsrMagTag = tagNames.waletensormag;
        if( !factory.have_entry( waleTsrMagTag ) ){
          typedef WaleTensorMagnitude::Builder waleStrTsrMagT;
          factory.register_expression( scinew waleStrTsrMagT(waleTsrMagTag, velTags ) );
        }
      }
        break;
        
        // ---------------------------------------------------------------------
      case DYNAMIC: {
        strTsrMagTag = tagNames.straintensormag;//( "StrainTensorMagnitude", Expr::STATE_NONE );

        Expr::TagList dynamicSmagTagList;
        dynamicSmagTagList.push_back( strTsrMagTag );
        dynamicSmagTagList.push_back( tagNames.dynamicsmagcoef);

        // if the DYNAMIC model is turned on, then create an expression for the dynamic smagorinsky coefficient
        dynSmagCoefTag = tagNames.dynamicsmagcoef;
        
        if( !factory.have_entry( dynSmagCoefTag )&&
            !factory.have_entry( strTsrMagTag )     ){
          typedef DynamicSmagorinskyCoefficient::Builder dynSmagConstT;
          factory.register_expression( scinew dynSmagConstT(dynamicSmagTagList,
                                                            velTags,
                                                            densTag,
                                                            isConstDensity) );
        }
        
      }
        break;

        // ---------------------------------------------------------------------
      default:
        break;
    }

    if( !factory.have_entry( turbViscTag ) ){
      // NOTE: You may need to cleave the turbulent viscosity from its parents
      // in case you run into problems with your simulation. The default behavior
      // of Wasatch is to extrapolate the turbulent viscosity at all patch boundaries.
      // If this extrapolation leads to problems you should consider excluding
      // the extrapolation and communicating the TurbulentViscosity instead.
      // To get rid of extrapolation, go to the TurbulentViscosity.cc expression
      // and comment out the "exOp_->apply_to_field(result)" line at the end
      // of the evaluate method.
      typedef TurbulentViscosity::Builder TurbViscT;
      factory.register_expression( scinew TurbViscT(turbViscTag, densTag, strTsrMagTag, waleTsrMagTag, vremanTsrMagTag, dynSmagCoefTag, turbParams ) );
//      const Expr::ExpressionID turbViscID = factory.register_expression( scinew TurbViscT(turbViscTag, densTag, strTsrMagTag, waleTsrMagTag, vremanTsrMagTag, dynSmagCoefTag, turbParams ) );
//      factory.cleave_from_parents(turbViscID);
    }
  }
  
  //==================================================================

  // note that the ordering of Vel1T and Vel2T are very important, and
  // must be consistent with the order of the velocity tags passed
  // into the Strain constructor.
  template< typename FaceT > struct StrainHelper;
  // nomenclature: XSurfXField - first letter is volume type: S, X, Y, Z
  // then it is followed by the field type
  template<> struct StrainHelper<SpatialOps::structured::XSurfXField>
  {
    // XSurfXField - XVol-XSurf
    // tau_xx
    typedef XVolField Vel1T;
    typedef XVolField Vel2T;
  };
  template<> struct StrainHelper<SpatialOps::structured::XSurfYField>
  {
    // XSurfYField - XVol-YSurf
    // tau_yx (tau on a y face in the x direction)
    typedef XVolField Vel1T;
    typedef YVolField Vel2T;
  };
  template<> struct StrainHelper<SpatialOps::structured::XSurfZField>
  {
    // XSurfZField - XVol-ZSurf
    // tau_zx (tau on a z face in the x direction)
    typedef XVolField Vel1T;
    typedef ZVolField Vel2T;
  };

  template<> struct StrainHelper<SpatialOps::structured::YSurfXField>
  {
    // tau_xy
    typedef YVolField Vel1T;
    typedef XVolField Vel2T;
  };
  template<> struct StrainHelper<SpatialOps::structured::YSurfYField>
  {
    // tau_yy
    typedef YVolField Vel1T;
    typedef YVolField Vel2T;
  };
  template<> struct StrainHelper<SpatialOps::structured::YSurfZField>
  {
    // tau_zy
    typedef YVolField Vel1T;
    typedef ZVolField Vel2T;
  };

  template<> struct StrainHelper<SpatialOps::structured::ZSurfXField>
  {
    // tau_xz
    typedef ZVolField Vel1T;
    typedef XVolField Vel2T;
  };
  template<> struct StrainHelper<SpatialOps::structured::ZSurfYField>
  {
    // tau_yz
    typedef ZVolField Vel1T;
    typedef YVolField Vel2T;
  };
  template<> struct StrainHelper<SpatialOps::structured::ZSurfZField>
  {
    // tau_zz
    typedef ZVolField Vel1T;
    typedef ZVolField Vel2T;
  };

  //==================================================================

  template< typename FieldT> struct NormalFaceSelector;

  template<> struct NormalFaceSelector<SpatialOps::structured::XVolField>
  {
  private:
    typedef SpatialOps::structured::XVolField FieldT;
  public:
    typedef SpatialOps::structured::FaceTypes<FieldT>::XFace NormalFace;
  };

  template<> struct NormalFaceSelector<SpatialOps::structured::YVolField>
  {
  private:
    typedef SpatialOps::structured::YVolField FieldT;
  public:
    typedef SpatialOps::structured::FaceTypes<FieldT>::YFace NormalFace;
  };

  template<> struct NormalFaceSelector<SpatialOps::structured::ZVolField>
  {
  private:
    typedef SpatialOps::structured::ZVolField FieldT;
  public:
    typedef SpatialOps::structured::FaceTypes<FieldT>::ZFace NormalFace;
  };

  //==================================================================

  Expr::Tag mom_tag( const std::string& momName )
  {
    return Expr::Tag( momName, Expr::STATE_N );
  }

  //==================================================================

  Expr::Tag rhs_part_tag( const Expr::Tag& momTag )
  {
    return Expr::Tag( momTag.name() + "_rhs_partial", Expr::STATE_NONE );
  }

  
  //==================================================================
  
  void set_vel_star_tags( Expr::TagList velTags,
                         Expr::TagList& velStarTags )
  {
    const TagNames& tagNames = TagNames::self();
    if( velTags[0] != Expr::Tag() ) velStarTags.push_back( Expr::Tag(velTags[0].name() + tagNames.star, Expr::STATE_NONE) );
    else         velStarTags.push_back( Expr::Tag() );
    if( velTags[1] != Expr::Tag() ) velStarTags.push_back( Expr::Tag(velTags[1].name() + tagNames.star, Expr::STATE_NONE) );
    else         velStarTags.push_back( Expr::Tag() );
    if( velTags[2] != Expr::Tag() ) velStarTags.push_back( Expr::Tag(velTags[2].name() + tagNames.star, Expr::STATE_NONE) );
    else         velStarTags.push_back( Expr::Tag() );
  }
  
  //==================================================================
  
  void set_vel_tags( Uintah::ProblemSpecP params,
                    Expr::TagList& velTags )
  {
    std::string xvelname, yvelname, zvelname;
    Uintah::ProblemSpecP doxvel,doyvel,dozvel;
    doxvel = params->get( "X-Velocity", xvelname );
    doyvel = params->get( "Y-Velocity", yvelname );
    dozvel = params->get( "Z-Velocity", zvelname );
    if( doxvel ) velTags.push_back( Expr::Tag(xvelname, Expr::STATE_NONE) );
    else         velTags.push_back( Expr::Tag() );
    if( doyvel ) velTags.push_back( Expr::Tag(yvelname, Expr::STATE_NONE) );
    else         velTags.push_back( Expr::Tag() );
    if( dozvel ) velTags.push_back( Expr::Tag(zvelname, Expr::STATE_NONE) );
    else         velTags.push_back( Expr::Tag() );
  }
  
  //==================================================================
  
  void set_mom_tags( Uintah::ProblemSpecP params,
                    Expr::TagList& momTags )
  {
    std::string xmomname, ymomname, zmomname;
    Uintah::ProblemSpecP doxmom,doymom,dozmom;
    doxmom = params->get( "X-Momentum", xmomname );
    doymom = params->get( "Y-Momentum", ymomname );
    dozmom = params->get( "Z-Momentum", zmomname );
    if( doxmom ) momTags.push_back( Expr::Tag(xmomname, Expr::STATE_N) );
    else         momTags.push_back( Expr::Tag() );
    if( doymom ) momTags.push_back( Expr::Tag(ymomname, Expr::STATE_N) );
    else         momTags.push_back( Expr::Tag() );
    if( dozmom ) momTags.push_back( Expr::Tag(zmomname, Expr::STATE_N) );
    else         momTags.push_back( Expr::Tag() );
  }
  
  //==================================================================
  
  template< typename FieldT >
  void
  set_tau_tags( const bool* doMom,
               const bool isViscous,
               Expr::TagList& tauTags )
  {
    const Direction stagLoc = get_staggered_location<FieldT>();
    std::string thisMomDirName;
    switch (stagLoc) {
      case XDIR:
        thisMomDirName = "x";
        break;
      case YDIR:
        thisMomDirName = "y";
        break;
      case ZDIR:
        thisMomDirName = "z";
        break;
      case NODIR:
      default:
        thisMomDirName = "";
        break;
    }

    if( doMom[0] && isViscous ) tauTags.push_back( Expr::Tag("tau_x" + thisMomDirName , Expr::STATE_NONE) );
    else                        tauTags.push_back( Expr::Tag() );
    if( doMom[1] && isViscous ) tauTags.push_back( Expr::Tag("tau_y" + thisMomDirName , Expr::STATE_NONE) );
    else                        tauTags.push_back( Expr::Tag() );
    if( doMom[2] && isViscous ) tauTags.push_back( Expr::Tag("tau_z" + thisMomDirName , Expr::STATE_NONE) );
    else                        tauTags.push_back( Expr::Tag() );
  }
  
  //==================================================================
  
  void set_convflux_tags( const bool* doMom,
                         Expr::TagList& cfTags,
                         const Expr::Tag thisMomTag )
  {
    if( doMom[0] ) cfTags.push_back( Expr::Tag(thisMomTag.name() + "_convFlux_x", Expr::STATE_NONE) );
    else         cfTags.push_back( Expr::Tag() );
    if( doMom[1] ) cfTags.push_back( Expr::Tag(thisMomTag.name() + "_convFlux_y", Expr::STATE_NONE) );
    else         cfTags.push_back( Expr::Tag() );
    if( doMom[2] ) cfTags.push_back( Expr::Tag(thisMomTag.name() + "_convFlux_z", Expr::STATE_NONE) );
    else         cfTags.push_back( Expr::Tag() );
  }

  //==================================================================

  /**
   *  \brief Register the Strain expression for the given face field
   */
  template< typename FaceFieldT >
  Expr::ExpressionID
  setup_strain( const Expr::Tag& strainTag,
                const Expr::Tag& vel1Tag,
                const Expr::Tag& vel2Tag,
                const Expr::Tag& dilTag,
                Expr::ExpressionFactory& factory )
  {
    typedef typename StrainHelper<FaceFieldT>::Vel1T Vel1T;  // type of velocity component 1
    typedef typename StrainHelper<FaceFieldT>::Vel2T Vel2T;  // type of velocity component 2
    typedef SVolField                                ViscT;  // type of viscosity

    typedef typename Strain< FaceFieldT, Vel1T, Vel2T >::Builder StrainT;

    return factory.register_expression( scinew StrainT( strainTag, vel1Tag, vel2Tag, dilTag ) );
  }

  //==================================================================

  template< typename FieldT >
  Expr::ExpressionID
  register_strain_tensor( const bool* const doMom,
                         const bool isViscous,
                         const Expr::TagList& velTags,
                         Expr::TagList& tauTags,
                         const Expr::Tag& dilTag,
                         Expr::ExpressionFactory& factory )
  {
    const Direction stagLoc = get_staggered_location<FieldT>();

    typedef typename SpatialOps::structured::FaceTypes<FieldT>::XFace XFace;
    typedef typename SpatialOps::structured::FaceTypes<FieldT>::YFace YFace;
    typedef typename SpatialOps::structured::FaceTypes<FieldT>::ZFace ZFace;

    set_tau_tags<FieldT>( doMom, isViscous, tauTags );
    const Expr::Tag& tauxt = tauTags[0];
    const Expr::Tag& tauyt = tauTags[1];
    const Expr::Tag& tauzt = tauTags[2];
    
    Expr::ExpressionID normalStrainID;
    
    const int thisVelIdx = (stagLoc == XDIR) ? 0 : ( (stagLoc == YDIR) ? 1 : 2 );    
    const Expr::Tag& thisVelTag = velTags[thisVelIdx];

    // register necessary strain expression when the flow is viscous
    if ( isViscous ) {
      if( doMom[0] ){
        const Expr::ExpressionID strainID = setup_strain< XFace >( tauxt, thisVelTag, velTags[0], dilTag, factory );
        if( stagLoc == XDIR )  normalStrainID = strainID;
      }
      if( doMom[1] ){
        const Expr::ExpressionID strainID = setup_strain< YFace >( tauyt, thisVelTag, velTags[1], dilTag, factory );
        if( stagLoc == YDIR )  normalStrainID = strainID;
      }
      if( doMom[2] ){
        const Expr::ExpressionID strainID = setup_strain< ZFace >( tauzt, thisVelTag, velTags[2], dilTag, factory );
        if( stagLoc == ZDIR )  normalStrainID = strainID;
      }
      factory.cleave_from_children( normalStrainID );
      factory.cleave_from_parents( normalStrainID  );
    }
    return normalStrainID;
  }
  
  //==================================================================

  template< typename FluxT, typename AdvelT >
  Expr::ExpressionID
  setup_convective_flux( const Expr::Tag& fluxTag,
                         const Expr::Tag& momTag,
                         const Expr::Tag& advelTag, Expr::ExpressionFactory& factory )
  {
    typedef typename SpatialOps::structured::VolType<FluxT>::VolField  MomT;
    typedef typename SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Interpolant, MomT,   FluxT >::type  MomInterpOp;
    typedef typename SpatialOps::structured::OperatorTypeBuilder< SpatialOps::Interpolant, AdvelT, FluxT >::type  AdvelInterpOp;
    typedef typename ConvectiveFlux<MomInterpOp, AdvelInterpOp >::Builder ConvFlux;
    return factory.register_expression( scinew ConvFlux( fluxTag, momTag, advelTag ) );
  }

  //==================================================================
  
  template< typename FieldT >
  Expr::ExpressionID
  register_convective_fluxes( const bool* const doMom,
                              const Expr::TagList& velTags,
                              Expr::TagList& cfTags,
                              const Expr::Tag& momTag,
                              Expr::ExpressionFactory& factory )
  {
    set_convflux_tags( doMom, cfTags, momTag );
    const Expr::Tag cfxt = cfTags[0];
    const Expr::Tag cfyt = cfTags[1];
    const Expr::Tag cfzt = cfTags[2];

    typedef typename SpatialOps::structured::FaceTypes<FieldT>::XFace XFace;
    typedef typename SpatialOps::structured::FaceTypes<FieldT>::YFace YFace;
    typedef typename SpatialOps::structured::FaceTypes<FieldT>::ZFace ZFace;

    Expr::ExpressionID normalConvFluxID;
    Direction stagLoc = get_staggered_location<FieldT>();
    
    if( doMom[0] ){
      const Expr::ExpressionID id = setup_convective_flux< XFace, XVolField >( cfxt, momTag, velTags[0], factory );
      if( stagLoc == XDIR )  normalConvFluxID = id;
    }
    if( doMom[1] ){
      const Expr::ExpressionID id = setup_convective_flux< YFace, YVolField >( cfyt, momTag, velTags[1], factory );
      if( stagLoc == YDIR )  normalConvFluxID = id;
    }
    if( doMom[2] ){
      const Expr::ExpressionID id = setup_convective_flux< ZFace, ZVolField >( cfzt, momTag, velTags[2], factory );
      if( stagLoc == ZDIR )  normalConvFluxID = id;
    }
    // convective fluxes require ghost updates after they are calculated
    // jcs note that we need to set BCs on these quantities as well.
    factory.cleave_from_children( normalConvFluxID );
    factory.cleave_from_parents ( normalConvFluxID );
    return normalConvFluxID;
  }

  //==================================================================

  template< typename FieldT >
  Expr::ExpressionID
  MomentumTransportEquation<FieldT>::
  get_mom_rhs_id( Expr::ExpressionFactory& factory,
                  const std::string velName,
                  const std::string momName,
                  Uintah::ProblemSpecP params,
                  const bool hasEmbeddedGeometry,
                  Uintah::SolverInterface& linSolver )
  {
    const Expr::Tag momTag = mom_tag( momName );
    const Expr::Tag rhsFullTag( momTag.name() + "_rhs_full", Expr::STATE_NONE );
    bool  enablePressureSolve = !(params->findBlock("DisablePressureSolve"));
    
    Expr::Tag volFracTag = Expr::Tag();
    Direction stagLoc = get_staggered_location<FieldT>();
    VolFractionNames& vNames = VolFractionNames::self();
    switch (stagLoc) {
      case XDIR:
        if (hasEmbeddedGeometry) volFracTag = vNames.xvol_frac_tag();
        break;
      case YDIR:
        if (hasEmbeddedGeometry) volFracTag = vNames.yvol_frac_tag();
        break;
      case ZDIR:
        if (hasEmbeddedGeometry) volFracTag = vNames.zvol_frac_tag();
        break;
      default:
        break;
    }    
    return factory.register_expression( new typename MomRHS<FieldT>::Builder( rhsFullTag, (enablePressureSolve ? pressure_tag() : Expr::Tag()), rhs_part_tag(momTag) , volFracTag ) );
  }

  //==================================================================

  template< typename FieldT >
  MomentumTransportEquation<FieldT>::
  MomentumTransportEquation( const std::string velName,
                             const std::string momName,
                             const Expr::Tag densTag,
                             const bool isConstDensity,
                             const Expr::Tag bodyForceTag,
                             const Expr::Tag srcTermTag,
                             GraphHelper& graphHelper,
                             Uintah::ProblemSpecP params,
                             TurbulenceParameters turbulenceParams,
                             const bool hasEmbeddedGeometry,
                             const bool hasMovingGeometry,
                             const Expr::ExpressionID rhsID,
                             Uintah::SolverInterface& linSolver,
                             Uintah::SimulationStateP sharedState)
    : TransportEquation( momName,
                         rhsID,
                         get_staggered_location<FieldT>(),
                         isConstDensity,
                         hasEmbeddedGeometry,
                         params ),
      isViscous_       ( params->findBlock("Viscosity") ? true : false ),
      isConstDensity_  ( isConstDensity                       ),
      isTurbulent_     ( turbulenceParams.turbModelName != NOTURBULENCE ),
      thisVelTag_      ( Expr::Tag(velName, Expr::STATE_NONE) ),
      densityTag_      ( densTag                              ),
      normalStrainID_  ( Expr::ExpressionID::null_id()        ),
      normalConvFluxID_( Expr::ExpressionID::null_id()        ),
      pressureID_      ( Expr::ExpressionID::null_id()        )
  {
    solverParams_ = NULL;
    set_vel_tags( params, velTags_ );

    Expr::ExpressionFactory& factory = *(graphHelper.exprFactory);

    thisMomName_ = momName;
    const Expr::Tag thisMomTag = mom_tag( thisMomName_ );
    
    const TagNames& tagNames = TagNames::self();
    
    std::string xmomname, ymomname, zmomname; // these are needed to construct fx, fy, and fz for pressure RHS
    bool doMom[3];
    doMom[0] = params->get( "X-Momentum", xmomname );
    doMom[1] = params->get( "Y-Momentum", ymomname );
    doMom[2] = params->get( "Z-Momentum", zmomname );

    //_____________
    // volume fractions for embedded boundaries Terms
    VolFractionNames& vNames = VolFractionNames::self();
    Expr::Tag volFracTag =     this->has_embedded_geometry()               ? vNames.svol_frac_tag() : Expr::Tag();
    Expr::Tag xAreaFracTag = ( this->has_embedded_geometry() && doMom[0] ) ? vNames.xvol_frac_tag() : Expr::Tag();
    Expr::Tag yAreaFracTag = ( this->has_embedded_geometry() && doMom[1] ) ? vNames.yvol_frac_tag() : Expr::Tag();
    Expr::Tag zAreaFracTag = ( this->has_embedded_geometry() && doMom[2] ) ? vNames.zvol_frac_tag() : Expr::Tag();

    //__________________
    // convective fluxes
    Expr::TagList cfTags; // these tags will be filled by register_convective_fluxes
    normalConvFluxID_ = register_convective_fluxes<FieldT>(doMom, velTags_, cfTags, thisMomTag, factory);

    //__________________
    // dilatation - needed by pressure source term and strain tensor
    const Expr::Tag dilTag = tagNames.dilatation;
    if( !factory.have_entry( dilTag ) ){
      typedef typename Dilatation<SVolField,XVolField,YVolField,ZVolField>::Builder Dilatation;
      // if dilatation expression has not been registered, then register it
      factory.register_expression( new Dilatation(dilTag, velTags_) );
    }

    //___________________________________
    // diffusive flux (strain components)
    Expr::TagList tauTags;
    normalStrainID_ = register_strain_tensor<FieldT>(doMom, isViscous_, velTags_, tauTags, dilTag, factory);
    
    //--------------------------------------
    // TURBULENCE
    // check if we have a turbulence model turned on
    // check if the flow is viscous
    const Expr::Tag viscTag = (isViscous_) ? parse_nametag( params->findBlock("Viscosity")->findBlock("NameTag") ) : Expr::Tag();
    
    bool enableTurbulenceModel = !(params->findBlock("DisableTurbulenceModel"));
    const Expr::Tag turbViscTag = tagNames.turbulentviscosity;
    if ( isTurbulent_ && isViscous_ && enableTurbulenceModel ) {
      register_turbulence_expressions(turbulenceParams, factory, velTags_, densTag, is_constant_density() );
      factory.attach_dependency_to_expression(turbViscTag, viscTag);
    }
    // END TURBULENCE
    //--------------------------------------

    //_________________________________________________________
    // partial rhs:
    // register expression to calculate the partial RHS (absent
    // pressure gradient) for use in the projection
    Expr::Tag volTag = Expr::Tag();
    switch (stagLoc_) {
      case XDIR:
        volTag = xAreaFracTag;
        break;
      case YDIR:
        volTag = yAreaFracTag;
        break;
      case ZDIR:
        volTag = zAreaFracTag;
        break;        
      default:
        break;
    }
    const Expr::ExpressionID momRHSPartID = factory.register_expression(
        new typename MomRHSPart<FieldT>::Builder( rhs_part_tag( thisMomTag ),
                                                  cfTags[0] , cfTags[1] , cfTags[2] , viscTag,
                                                  tauTags[0], tauTags[1], tauTags[2], densityTag_,
                                                  bodyForceTag, srcTermTag,
                                                  volTag) );
    factory.cleave_from_parents ( momRHSPartID );
    
    //__________________
    // drhou/dt - needed by pressure BCs
    // Here we calculate drhoudt using rhou and rhou_old. The catch is to remember
    // to set the momentum_rhs_full to zero at wall/inlet boundaries.
    // One could also calculate this time derivative from rho, u and rho_old, u_old.
    // We may consider this option later.
    bool enabledudtInPRHS = !(params->findBlock("Disabledmomdt"));    
    if (enabledudtInPRHS) {
      OldVariable& oldVar = OldVariable::self();
      Expr::Tag dthisMomdtTag = Expr::Tag( "d_" + thisMomName_ + "_dt" , Expr::STATE_NONE );
      Expr::Tag thisMomOldTag = Expr::Tag( thisMomName_  + "_old", Expr::STATE_NONE );
      Expr::Tag thisMomOldOldTag = Expr::Tag( thisMomName_  + "_old_old", Expr::STATE_NONE );
      oldVar.add_variable<FieldT>( ADVANCE_SOLUTION, thisMomTag);
      oldVar.add_variable<FieldT>( ADVANCE_SOLUTION, thisMomOldTag);
      factory.register_expression( new typename TimeDerivative<FieldT>::Builder(dthisMomdtTag,thisMomOldTag,thisMomOldOldTag,tagNames.timestep));
    }

    //__________________
    // Pressure source term
    if (!isConstDensity) {
      // calculating velocity at the next time step    
      Expr::Tag thisVelStarTag = Expr::Tag( thisVelTag_.name() + tagNames.star, Expr::STATE_NONE);
      Expr::Tag convTermWeak   = Expr::Tag( thisVelTag_.name() + "_weak_convective_term", Expr::STATE_NONE);
      if( !factory.have_entry( thisVelStarTag ) ){
        OldVariable& oldPressure = OldVariable::self();
        oldPressure.add_variable<SVolField>( ADVANCE_SOLUTION, pressure_tag() );
        const Expr::Tag oldPressureTag = Expr::Tag (pressure_tag().name() + "_old", Expr::STATE_NONE);
        convTermWeakID_ = factory.register_expression( new typename WeakConvectiveTerm<FieldT>::Builder( convTermWeak, thisVelTag_, velTags_));
//        factory.cleave_from_children( convTermWeakID_ );
        factory.cleave_from_parents ( convTermWeakID_ );
        factory.register_expression( new typename VelEst<FieldT>::Builder( thisVelStarTag, thisVelTag_, convTermWeak, tauTags, densTag, viscTag, oldPressureTag, tagNames.timestep ));
      }
    }
    
    Expr::Tag pSourceTag = Expr::Tag( "pressure-source-term", Expr::STATE_NONE);
    if( !factory.have_entry( pSourceTag ) ){
      Expr::Tag densStarTag = Expr::Tag(densTag.name() + tagNames.star, Expr::CARRY_FORWARD);
      Expr::Tag dens2StarTag = Expr::Tag(densTag.name() + tagNames.doubleStar, Expr::CARRY_FORWARD);
      Expr::TagList velStarTags = Expr::TagList();
      
      set_vel_star_tags( velTags_, velStarTags );
      set_mom_tags( params, momTags_ );
      
      // registering the expressiong for pressure source term
      factory.register_expression( new typename PressureSource::Builder( pSourceTag, momTags_, velStarTags, isConstDensity, densTag, densStarTag, dens2StarTag, dilTag, tagNames.timestep));
    }
    
    
    //__________________
    // calculate velocity at the current time step    
    factory.register_expression( new typename PrimVar<FieldT,SVolField>::Builder( thisVelTag_, thisMomTag, densityTag_, volTag ));
    
    //__________________
    // pressure
    bool enablePressureSolve = !(params->findBlock("DisablePressureSolve"));

    if (enablePressureSolve) {
      if( !factory.have_entry( pressure_tag() ) ){
        Uintah::ProblemSpecP pressureParams = params->findBlock( "Pressure" );
        
        bool usePressureRefPoint = false;
        double refPressureValue = 0.0;
        SCIRun::IntVector refPressureLocation(0,0,0);
        if (pressureParams->findBlock("ReferencePressure")) {
          usePressureRefPoint = true;
          Uintah::ProblemSpecP refPressureParams = pressureParams->findBlock("ReferencePressure");
          refPressureParams->getAttribute("value", refPressureValue);
          refPressureParams->get("ReferenceCell", refPressureLocation);
        }
        
        bool use3DLaplacian = true;
        pressureParams->getWithDefault("Use3DLaplacian",use3DLaplacian, true);
        
        solverParams_ = linSolver.readParameters( pressureParams, "",
                                                 sharedState );
        solverParams_->setSolveOnExtraCells( false );
        solverParams_->setUseStencil4( true );
        solverParams_->setOutputFileName( "WASATCH" );
        
        // if pressure expression has not be registered, then register it
        Expr::Tag fxt, fyt, fzt;
        if( doMom[0] )  fxt = Expr::Tag( xmomname + "_rhs_partial", Expr::STATE_NONE );
        if( doMom[1] )  fyt = Expr::Tag( ymomname + "_rhs_partial", Expr::STATE_NONE );
        if( doMom[2] )  fzt = Expr::Tag( zmomname + "_rhs_partial", Expr::STATE_NONE );

        // add drhoudt term to the pressure rhs
        Expr::Tag dxmomdtt, dymomdtt, dzmomdtt;
        if (enabledudtInPRHS) {
          if( doMom[0] )  dxmomdtt = Expr::Tag( "d_" + xmomname + "_dt", Expr::STATE_NONE );
          if( doMom[1] )  dymomdtt = Expr::Tag( "d_" + ymomname + "_dt", Expr::STATE_NONE );
          if( doMom[2] )  dzmomdtt = Expr::Tag( "d_" + zmomname + "_dt", Expr::STATE_NONE );
        }

        Expr::TagList ptags;
        ptags.push_back( pressure_tag() );
        ptags.push_back( Expr::Tag( pressure_tag().name() + "_rhs", pressure_tag().context() ) );
        const Expr::ExpressionBuilder* const pbuilder = new typename Pressure::Builder( ptags, fxt, fyt, fzt, dxmomdtt, dymomdtt, dzmomdtt,
                                                                                        pSourceTag, tagNames.timestep, volFracTag, 
                                                                                        hasMovingGeometry, usePressureRefPoint, refPressureValue, 
                                                                                        refPressureLocation, use3DLaplacian,
                                                                                        *solverParams_, linSolver);
        pressureID_ = factory.register_expression( pbuilder );
        factory.cleave_from_children( pressureID_ );
        factory.cleave_from_parents ( pressureID_ );
      }
      else {
        pressureID_ = factory.get_id( pressure_tag() );
      }
    }
    
    // Kinetic energy calculation, if necessary
    if ( params->findBlock("CalculateKE") ) {
      Uintah::ProblemSpecP keSpec = params->findBlock("CalculateKE");
      bool isTotalKE = true;
      keSpec->getAttribute("total", isTotalKE);
      if (isTotalKE) { // calculate total kinetic energy. then follow that with a reduction variable
        if (!factory.have_entry( TagNames::self().totalKineticEnergy )) {
          bool outputKE = true;
          keSpec->getAttribute("output", outputKE);
          
          // we need to create two expressions
          const Expr::Tag tkeTempTag("TotalKE_temp", Expr::STATE_NONE);
          factory.register_expression(scinew typename TotalKineticEnergy<XVolField,YVolField,ZVolField>::Builder( tkeTempTag,
                                                                                                                 velTags_[0],velTags_[1],velTags_[2] ),true);
          
          ReductionHelper::self().add_variable<double, ReductionSumOpT>(ADVANCE_SOLUTION, TagNames::self().totalKineticEnergy, tkeTempTag, outputKE, false);
        }
      } else if (!factory.have_entry( TagNames::self().kineticEnergy )) { // calculate local, pointwise kinetic energy
        const Expr::ExpressionID keID = factory.register_expression(scinew typename KineticEnergy<SVolField,XVolField,YVolField,ZVolField>::Builder( TagNames::self().kineticEnergy,
                                                                                                                                                                 velTags_[0],velTags_[1],velTags_[2] ), true);
        graphHelper.rootIDs.insert( keID );
      }
    }

  }

  //------------------------------------------------------------------

  template< typename FieldT >
  MomentumTransportEquation<FieldT>::
  ~MomentumTransportEquation()
  {
    delete solverParams_;
  }

  //------------------------------------------------------------------

  template< typename FieldT >
  void
  MomentumTransportEquation<FieldT>::
  setup_initial_boundary_conditions( const GraphHelper& graphHelper,
                                     const Uintah::PatchSet* const localPatches,
                                     const PatchInfoMap& patchInfoMap,
                                     const Uintah::MaterialSubset* const materials,
                                     const std::map<std::string, std::set<std::string> >& bcFunctorMap)
  {
    Expr::ExpressionFactory& factory = *graphHelper.exprFactory;

    // multiply the initial condition by the volume fraction for embedded geometries
    if (hasEmbeddedGeometry_) {
      VolFractionNames& vNames = VolFractionNames::self();
      Expr::Tag volFracTag = Expr::Tag();
      switch (stagLoc_) {
        case XDIR:
          if (hasEmbeddedGeometry_) volFracTag = vNames.xvol_frac_tag();
          break;
        case YDIR:
          if (hasEmbeddedGeometry_) volFracTag = vNames.yvol_frac_tag();
          break;
        case ZDIR:
          if (hasEmbeddedGeometry_) volFracTag = vNames.zvol_frac_tag();
          break;
        default:
          break;
      }
      
      //create modifier expression
      typedef ExprAlgebra<FieldT> ExprAlgbr;
      Expr::TagList theTagList;
      theTagList.push_back(volFracTag);
      //theTagList.push_back(mom_tag(thisMomName_));
      Expr::Tag modifierTag = Expr::Tag( this->solution_variable_name() + "_init_cond_modifier", Expr::STATE_NONE);
      factory.register_expression( new typename ExprAlgbr::Builder(modifierTag,
                                                                   theTagList,
                                                                   ExprAlgbr::PRODUCT,
                                                                   true) );
      
      for( int ip=0; ip<localPatches->size(); ++ip ){
        
        // get the patch subset
        const Uintah::PatchSubset* const patches = localPatches->getSubset(ip);
        
        // loop over every patch in the patch subset
        for( int ipss=0; ipss<patches->size(); ++ipss ){
          
          // loop over materials
          for( int im=0; im<materials->size(); ++im ){
            //    if (hasVolFrac_) {
            // attach the modifier expression to the target expression
            factory.attach_modifier_expression( modifierTag, mom_tag(thisMomName_), true );
            //    }
            
          }
        }
      }
    }
        
    typedef typename NormalFaceSelector<FieldT>::NormalFace NormalFace;

    // set initial bcs for momentum
    if (factory.have_entry(mom_tag(thisMomName_))) {
      process_boundary_conditions<FieldT>( Expr::Tag( this->solution_variable_name(),
                                                      Expr::STATE_N ),
                                           this->solution_variable_name(),
                                           this->staggered_location(),
                                           graphHelper,
                                           localPatches,
                                           patchInfoMap,
                                           materials, bcFunctorMap );
    }

    // set bcs for velocity - cos we don't have a mechanism now to set them
    // on interpolated density field
    Expr::Tag velTag;
    switch (this->staggered_location()) {
      case XDIR:  velTag=velTags_[0];  break;
      case YDIR:  velTag=velTags_[1];  break;
      case ZDIR:  velTag=velTags_[2];  break;
      default:                         break;
    }
    if (factory.have_entry(velTag)) {
//      process_boundary_conditions<FieldT>( velTag,
//                                           velTag.name(),
//                                           this->staggered_location(),
//                                           graphHelper,
//                                           localPatches,
//                                           patchInfoMap,
//                                           materials );
    }

    // set bcs for pressure
    // We cannot set pressure BCs here using Wasatch's BC techniques because
    // we need to set the BCs AFTER the pressure solve. We had to create
    // a uintah task for that. See Pressure.cc
    
    // set bcs for partial rhs
    if (factory.have_entry(rhs_part_tag(mom_tag(thisMomName_)))) {
      process_boundary_conditions<FieldT>( rhs_part_tag(mom_tag(thisMomName_)),
                                           rhs_part_tag(mom_tag(thisMomName_)).name(),
                                           this->staggered_location(),
                                           graphHelper,
                                           localPatches,
                                           patchInfoMap,
                                           materials, bcFunctorMap );
    }
    
    if (!isConstDensity_) {
      // set bcs for density
      const Expr::Tag densTag( densityTag_.name(), Expr::STATE_NONE );
      process_boundary_conditions<SVolField>( densTag,
                                              densTag.name(),
                                              NODIR,
                                              graphHelper,
                                              localPatches,
                                              patchInfoMap,
                                              materials, bcFunctorMap );

      // set bcs for density_*
      const TagNames& tagNames = TagNames::self();
      const Expr::Tag densStarTag( densityTag_.name()+tagNames.star, Expr::STATE_NONE );
      const Expr::Tag densStarBCTag( densStarTag.name()+"_bc",Expr::STATE_NONE);
      Expr::ExpressionFactory& factory = *graphHelper.exprFactory;
      if (!factory.have_entry(densStarBCTag)){
        factory.register_expression ( new typename BCCopier<SVolField>::Builder(densStarBCTag, densTag) );
      }  
      process_boundary_conditions<SVolField>( densStarTag,
                                              densStarTag.name(),
                                              NODIR,
                                              graphHelper,
                                              localPatches,
                                              patchInfoMap,
                                              materials, bcFunctorMap,
                                              densTag.name(), 0, "Dirichlet", densStarBCTag.name() );
      
      // set bcs for velocity - cos we don't have a mechanism now to set them
      // on interpolated density field
      Expr::Tag velTag;
      switch (this->staggered_location()) {
        case XDIR:  velTag=velTags_[0];  break;
        case YDIR:  velTag=velTags_[1];  break;
        case ZDIR:  velTag=velTags_[2];  break;
        default:                         break;
      }
      process_boundary_conditions<FieldT>( velTag,
                                          velTag.name(),
                                           this->staggered_location(),
                                           graphHelper,
                                           localPatches,
                                           patchInfoMap,
                                           materials, bcFunctorMap );
    
    }
  }

  //------------------------------------------------------------------

  template< typename FieldT >
  void
  MomentumTransportEquation<FieldT>::
  setup_boundary_conditions( const GraphHelper& graphHelper,
                             const Uintah::PatchSet* const localPatches,
                             const PatchInfoMap& patchInfoMap,
                             const Uintah::MaterialSubset* const materials,
                             const std::map<std::string, std::set<std::string> >& bcFunctorMap)
  {
    typedef typename NormalFaceSelector<FieldT>::NormalFace NormalFace;

    // set bcs for momentum
    process_boundary_conditions<FieldT>( Expr::Tag( this->solution_variable_name(),
                                                    Expr::STATE_N ),
                                         this->solution_variable_name(),
                                         this->staggered_location(),
                                         graphHelper,
                                         localPatches,
                                         patchInfoMap,
                                         materials, bcFunctorMap );

    // set bcs for velocity - cos we don't have a mechanism now to set them
    // on interpolated density field
    Expr::Tag velTag;
    switch (this->staggered_location()) {
      case XDIR:  velTag=velTags_[0];  break;
      case YDIR:  velTag=velTags_[1];  break;
      case ZDIR:  velTag=velTags_[2];  break;
      default:                         break;
    }
    process_boundary_conditions<FieldT>( velTag,
                                         velTag.name(),
                                         this->staggered_location(),
                                         graphHelper,
                                         localPatches,
                                         patchInfoMap,
                                         materials, bcFunctorMap );
    if (!isConstDensity_) {
      // set bcs for density
      const Expr::Tag densTag( densityTag_.name(), Expr::STATE_NONE );
      process_boundary_conditions<SVolField>( densTag,
                                             densTag.name(),
                                             NODIR,
                                             graphHelper,
                                             localPatches,
                                             patchInfoMap,
                                             materials, bcFunctorMap );
      
      // set bcs for density_*
      const TagNames& tagNames = TagNames::self();
      const Expr::Tag densStarTag( densityTag_.name()+tagNames.star, Expr::CARRY_FORWARD );
      const Expr::Tag densStarBCTag( densStarTag.name()+"_bc",Expr::STATE_NONE);
      Expr::ExpressionFactory& factory = *graphHelper.exprFactory;
      if (!factory.have_entry(densStarBCTag)){
        factory.register_expression ( new typename BCCopier<SVolField>::Builder(densStarBCTag, densityTag_) );
      }  
      process_boundary_conditions<SVolField>( densStarTag,
                                             densStarTag.name(),
                                             NODIR,
                                             graphHelper,
                                             localPatches,
                                             patchInfoMap,
                                             materials, bcFunctorMap,
                                             densityTag_.name(), 0, "Dirichlet", densStarBCTag.name());
    }
    // set bcs for pressure
//    process_boundary_conditions<SVolField>( pressure_tag(),
//                                            "pressure",
//                                            NODIR,
//                                            graphHelper,
//                                            localPatches,
//                                            patchInfoMap,
//                                            materials );
    // set bcs for partial rhs
    process_boundary_conditions<FieldT>( rhs_part_tag(mom_tag(thisMomName_)),
                                         rhs_part_tag(mom_tag(thisMomName_)).name(),
                                         this->staggered_location(),
                                         graphHelper,
                                         localPatches,
                                         patchInfoMap,
                                         materials, bcFunctorMap );
    // set bcs for partial full rhs
    process_boundary_conditions<FieldT>( Expr::Tag(thisMomName_ + "_rhs_full", Expr::STATE_NONE),
                                        thisMomName_ + "_rhs_full",
                                        this->staggered_location(),
                                        graphHelper,
                                        localPatches,
                                        patchInfoMap,
                                        materials,bcFunctorMap );

    // set bcs for normal strains
    Expr::ExpressionFactory& factory = *graphHelper.exprFactory;
    if(isViscous_) {
      Expr::Tag normalStrainTag = factory.get_label(normalStrainID_);
      process_boundary_conditions<NormalFace>( normalStrainTag,
                                  normalStrainTag.name(),
                NODIR,
                graphHelper,
                localPatches,
                patchInfoMap,
                materials, bcFunctorMap);
    }

    // set bcs for normal convective fluxes
    Expr::Tag normalConvFluxTag = factory.get_label(normalConvFluxID_);
    process_boundary_conditions<NormalFace>( normalConvFluxTag,
                                normalConvFluxTag.name(),
                                NODIR,
                                graphHelper,
                                localPatches,
                                patchInfoMap,
                                materials, bcFunctorMap);

  }

  //------------------------------------------------------------------

  template< typename FieldT >
  Expr::ExpressionID
  MomentumTransportEquation<FieldT>::
  initial_condition( Expr::ExpressionFactory& icFactory )
  {
    // register an initial condition for da pressure
    Expr::Tag ptag(pressure_tag().name(), Expr::STATE_N);    
    if( !icFactory.have_entry( ptag ) ) {
      icFactory.register_expression( new typename Expr::ConstantExpr<SVolField>::Builder(ptag, 0.0 ) );
    }


    if( icFactory.have_entry( thisVelTag_ ) ) {
      typedef typename InterpolateExpression<SVolField, FieldT>::Builder Builder;
      Expr::Tag interpolatedDensityTag(densityTag_.name() +"_interp_" + this->dir_name(), Expr::STATE_NONE);
      icFactory.register_expression(scinew Builder(interpolatedDensityTag, Expr::Tag(densityTag_.name(),Expr::STATE_NONE)));
      
      // register expression to calculate the momentum initial condition from the initial conditions on
      // velocity and density in the cases that we are initializing velocity in the input file
      typedef ExprAlgebra<FieldT> ExprAlgbr;
      Expr::TagList theTagList;
      theTagList.push_back(thisVelTag_);
      theTagList.push_back(interpolatedDensityTag);
      return icFactory.register_expression( new typename ExprAlgbr::Builder( mom_tag(thisMomName_),
                                                                             theTagList,
                                                                             ExprAlgbr::PRODUCT ) );
    }
    
    return icFactory.get_id( Expr::Tag( this->solution_variable_name(), Expr::STATE_N ) );
  }

  //------------------------------------------------------------------

  //==================================================================
  // Explicit template instantiation
  template class MomentumTransportEquation< XVolField >;
  template class MomentumTransportEquation< YVolField >;
  template class MomentumTransportEquation< ZVolField >;
  
#define INSTANTIATE_SETUP_STRAIN(VOLT) \
  template Expr::ExpressionID setup_strain< SpatialOps::structured::FaceTypes<VOLT>::XFace > ( const Expr::Tag& strainTag, \
                                                                                                const Expr::Tag& vel1Tag,\
                                                                                                const Expr::Tag& vel2Tag,\
                                                                                                const Expr::Tag& dilTag,\
                                                                                                Expr::ExpressionFactory& factory ); \
  template Expr::ExpressionID setup_strain< SpatialOps::structured::FaceTypes<VOLT>::YFace > ( const Expr::Tag& strainTag, \
                                                                                                const Expr::Tag& vel1Tag,\
                                                                                                const Expr::Tag& vel2Tag,\
                                                                                                const Expr::Tag& dilTag,\
                                                                                                Expr::ExpressionFactory& factory ); \
  template Expr::ExpressionID setup_strain< SpatialOps::structured::FaceTypes<VOLT>::ZFace > ( const Expr::Tag& strainTag, \
                                                                                                const Expr::Tag& vel1Tag,\
                                                                                                const Expr::Tag& vel2Tag,\
                                                                                                const Expr::Tag& dilTag,\
                                                                                                Expr::ExpressionFactory& factory );
  INSTANTIATE_SETUP_STRAIN(XVolField);
  INSTANTIATE_SETUP_STRAIN(YVolField);
  INSTANTIATE_SETUP_STRAIN(ZVolField);
  //==================================================================

} // namespace Wasatch
