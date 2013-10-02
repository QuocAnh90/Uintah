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

#include "FluxLimiterInterpolant.h"
#include "OperatorTypes.h"
#include <CCA/Components/Wasatch/FieldAdaptor.h>

#include <cmath>

#include "spatialops/SpatialOpsDefs.h"
#include "spatialops/structured/FVTools.h"

#define LIMITER_TO_MACRO(limiterType) \
switch (limiterType) { \
  case Wasatch::SUPERBEE: \
    SUPERBEE \
  break;\
}

#define UPWIND(r) \
0.0

#define SUPERBEE(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, max(0.0, max( min(2.0*r, 1.0), min(r, 2.0) ) ) ) \
      ( 2.0 )

#define CHARM(r) \
cond  ( r < 0.0, 0.0 ) \
      ( r < infinity_, r * (3.0*r + 1)/( (r + 1)*(r + 1) ) ) \
      ( 3.0 )

#define KOREN(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, max( 0.0, min(min(2.0*r, (r + 2.0)/3.0) ,2.0 ) ) )\
      ( 2.0 )

#define MC(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, max( 0.0, min( min(2.0*r, 0.5*(r + 1.0)) ,2.0 ) ) )\
      ( 2.0 )

#define OSPRE(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, 1.5*(r*r + r)/(r*r + r + 1.0) )\
      ( 1.5 )

#define SMART(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, max( 0.0, min( min(2.0*r, 0.75*r + 0.25), 4.0 ) ) )\
      ( 4.0 )

#define VANLEER(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, (r + abs(r))/(1 + abs(r)) )\
      ( 2.0 )

#define HCUS(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, ( 1.5*(r + abs(r)) )/(r + 2.0) )\
      ( 3.0 )

#define MINMOD(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, max( 0.0, min(1.0, r) ) )\
      ( 1.0 )

#define HQUICK(r) \
cond  ( r < -infinity_, 0.0 ) \
      ( r < infinity_, ( 2*(r + abs(r)) )/(r + 3.0) )\
      (4.0 )

// vel is advective velocity, rm is gradient on minus side, rp is gradient on plus side
#define CALCULATE_INTERIOR_LIMITER(vel, rm, rp, LIMITERTYPE)  \
cond  ( vel > 0.0, LIMITERTYPE(rm) ) \
      ( vel < 0.0, LIMITERTYPE(rp) ) \
      ( 1.0 );

#define CALCULATE_BOUNDARY_LIMITER(flowDir, r, LIMITERTYPE)  \
cond  ( flowDir > 0.0, LIMITERTYPE(r) )                     \
      ( flowDir < 0.0, cond ( isBoundaryFace, 1.0)       \
                            ( 0.0) )                     \
      ( 1.0 );


//--------------------------------------------------------------------

template< typename PhiVolT, typename PhiFaceT >
FluxLimiterInterpolant<PhiVolT,PhiFaceT>::
FluxLimiterInterpolant( const std::vector<int>& dim,
                        const std::vector<bool> hasPlusFace,
                        const std::vector<bool> hasMinusBoundary )
: hasPlusBoundary_ (false),
  hasMinusBoundary_(false)
{
  
  const size_t direction = PhiFaceT::Location::FaceDir::value;
  
  switch (direction) {
    case SpatialOps::XDIR::value:
      unitNormal_ = SpatialOps::structured::IntVec(1,0,0);
      hasMinusBoundary_ = hasMinusBoundary[0];
      hasPlusBoundary_ = hasPlusFace[0];
      break;
    case SpatialOps::YDIR::value:
      unitNormal_ = SpatialOps::structured::IntVec(0,1,0);
      hasMinusBoundary_ = hasMinusBoundary[1];
      hasPlusBoundary_ = hasPlusFace[1];
      break;
    case SpatialOps::ZDIR::value:
      unitNormal_ = SpatialOps::structured::IntVec(0,0,1);
      hasMinusBoundary_ = hasMinusBoundary[2];
      hasPlusBoundary_ = hasPlusFace[2];
      break;
    default:
      unitNormal_ = SpatialOps::structured::IntVec(0,0,0);
      break;
  }
  
}

//--------------------------------------------------------------------

template< typename PhiVolT, typename PhiFaceT >
void
FluxLimiterInterpolant<PhiVolT,PhiFaceT>::
set_advective_velocity( const PhiFaceT &theAdvectiveVelocity )
{
# ifdef ENABLE_THREADS
  /*
   * Because this operator may be accessed by multiple expressions simultaneously,
   * there is the possibility that the advective velocity could be set by multiple
   * threads simultaneously.  Therefore, we lock this until the apply_to_field()
   * is done with the advective velocity to prevent race conditions.
   */
  mutex_.lock();
#endif
  advectiveVelocity_ = &theAdvectiveVelocity;
}

//--------------------------------------------------------------------

template< typename PhiVolT, typename PhiFaceT >
void
FluxLimiterInterpolant<PhiVolT,PhiFaceT>::
set_flux_limiter_type( Wasatch::ConvInterpMethods limiterType )
{
  limiterType_ = limiterType;
}

//--------------------------------------------------------------------

template< typename PhiVolT, typename PhiFaceT >
FluxLimiterInterpolant<PhiVolT,PhiFaceT>::
~FluxLimiterInterpolant()
{}

//--------------------------------------------------------------------

template< typename PhiVolT, typename PhiFaceT >
void
FluxLimiterInterpolant<PhiVolT,PhiFaceT>::
apply_embedded_boundaries( const PhiVolT &src, PhiFaceT &dest ) const {
  
  using namespace SpatialOps;
  using namespace SpatialOps::structured;


  build_src_fields(src);
  
  const MemoryWindow& wdest = dest.window_with_ghost(); // used for velocity & interpolated phi
  IntVec destExtent = wdest.extent() - unitNormal_*3 - unitNormal_ * (hasPlusBoundary_ ? 1 : 0);
  IntVec destBaseOffset = wdest.offset() + unitNormal_*2;
  // this is the destination field value - always on the boundary
  const MemoryWindow wd( wdest.glob_dim(),
                         destBaseOffset,
                         destExtent,
                         wdest.has_bc(0), wdest.has_bc(1), wdest.has_bc(2) );

  const MemoryType dMemType = dest.memory_device_type();  // destination memory type
  const unsigned short int dDevIdx = dest.device_index(); // destination device index
  typename PhiFaceT::value_type* destVals = const_cast<typename PhiFaceT::value_type*>( dest.field_values(dMemType, dDevIdx) );

  const MemoryType advelMemType = advectiveVelocity_->memory_device_type();  // advel memory type
  const unsigned short int advelDevIdx = advectiveVelocity_->device_index(); // advel device index
  typename PhiFaceT::value_type* velVals = const_cast<typename PhiFaceT::value_type*>( advectiveVelocity_->field_values(advelMemType,advelDevIdx) );

  PhiVolT d( wd, destVals, ExternalStorage, dMemType, dDevIdx );
  const PhiVolT aVel( wd, velVals, ExternalStorage, advelMemType, advelDevIdx );

  const PhiVolT& vfracmm = *srcFields_[0];
  const PhiVolT& vfracpp = *srcFields_[3];
  
  d <<= cond( aVel > 0.0 && vfracmm == 0.0, 0.0 )
            ( aVel < 0.0 && vfracpp == 0.0, 0.0 )
            ( d );  
}

//--------------------------------------------------------------------

template< typename PhiVolT, typename PhiFaceT >
void
FluxLimiterInterpolant<PhiVolT,PhiFaceT>::
build_src_fields( const PhiVolT& src ) const
{
  using namespace SpatialOps;
  using namespace SpatialOps::structured;
  
  srcFields_.clear();

  // build source fields
  const MemoryWindow& wsrc = src.window_with_ghost();
  for (int i=0; i<4; i++) {
    // i = 0: minus-minus
    // i = 1: minus
    // i = 2: plus
    // i = 3: plus-plus
    //  below is a depiction of the nomenclature, cells and faces...
    //  | minus-minus |  minus  || plus  | plus-plus
    const MemoryWindow srcwin( wsrc.glob_dim(),
                               wsrc.offset() + unitNormal_*i,
                               wsrc.extent() - unitNormal_*3,
                               wsrc.has_bc(0), wsrc.has_bc(1), wsrc.has_bc(2) );
    
    const PhiVolT * field = new PhiVolT( srcwin, src );
    
    srcFields_.push_back( PhiVolTPtr(field) );
  }
}

//--------------------------------------------------------------------

template< typename PhiVolT, typename PhiFaceT >
void
FluxLimiterInterpolant<PhiVolT,PhiFaceT>::
apply_to_field( const PhiVolT &src, PhiFaceT &dest )
{
  // This will calculate the flux limiter function psi. The convective flux is
  // written as: phi_face = phi_lo - psi*(phi_lo - phi_hi) where
  // phi_lo is a low order interpolant (i.e. Upwind)
  // phi_hi is a high order interpolant (i.e. central).
  /* Algorithm: TSAAD - TODO - DESCRIBE ALGORITHM IN DETAIL
   * Loop over faces
   */
  using namespace SpatialOps;
  using namespace SpatialOps::structured;

  const MemoryWindow& wsrc  = src.window_with_ghost();
  const MemoryWindow& wdest = dest.window_with_ghost(); // used for velocity & interpolated phi
  
  const MemoryType dMemType = dest.memory_device_type();  // destination memory type
  const unsigned short int dDevIdx = dest.device_index(); // destination device index
  typename PhiFaceT::value_type* destVals = const_cast<typename PhiFaceT::value_type*>( dest.field_values(dMemType, dDevIdx) );

  const MemoryType advelMemType = advectiveVelocity_->memory_device_type();  // advel memory type
  const unsigned short int advelDevIdx = advectiveVelocity_->device_index(); // advel device index
  typename PhiFaceT::value_type* velVals = const_cast<typename PhiFaceT::value_type*>( advectiveVelocity_->field_values(advelMemType,advelDevIdx) );

  int pm[2]={1,-1}; // plus or minus face
  int zo[2]={0,1};  // zero and one
  const double infinity_ = 1.0e10;

  IntVec extent = wsrc.extent() - unitNormal_*wsrc.glob_dim() + unitNormal_;
  IntVec destExtent = wdest.extent() - unitNormal_*wdest.glob_dim() + unitNormal_;
  IntVec baseOffset;
  IntVec destBaseOffset;
  
  // start with patch boundaries
  for (int direc=0; direc<2; direc++) {
    baseOffset = wsrc.offset() + (unitNormal_*wsrc.glob_dim() - unitNormal_ )* zo[direc]; // src base offset
    destBaseOffset = wdest.offset() + (unitNormal_*wdest.glob_dim() - unitNormal_ - unitNormal_*hasPlusBoundary_ )* zo[direc] + unitNormal_*(1-zo[direc]); // destination base offset - depends on presence of plus boundary
    
    // this is the destination field value - always on the boundary
    const MemoryWindow wd( wdest.glob_dim(),
                           destBaseOffset,
                           destExtent,
                           wdest.has_bc(0), wdest.has_bc(1), wdest.has_bc(2) );
    
    // ghost cell: on a minus face, this is src-minus. on a plus face, this is src-plus
    const MemoryWindow ws1( wsrc.glob_dim(),
                            baseOffset,
                            extent,
                            wsrc.has_bc(0), wsrc.has_bc(1), wsrc.has_bc(2) );
    
    // first interior cell: on a minus face, this is src-plus. on a plus face this is src-minus
    const MemoryWindow ws2( wsrc.glob_dim(),
                            baseOffset + unitNormal_ * pm[direc],
                            extent,
                            wsrc.has_bc(0), wsrc.has_bc(1), wsrc.has_bc(2) );
        
    // second interior cell: on a minus face, this is src-plus-plus. on a plus face, this is src-minus-minus
    const MemoryWindow ws3( wsrc.glob_dim(),
                            baseOffset  + unitNormal_ * pm[direc] * 2,
                            extent,
                            wsrc.has_bc(0), wsrc.has_bc(1), wsrc.has_bc(2) );

    PhiVolT d( wd, destVals, ExternalStorage, dMemType, dDevIdx );
    const PhiVolT aVel( wd, velVals, ExternalStorage, advelMemType, advelDevIdx );
    const PhiVolT s1( ws1, src );
    const PhiVolT s2( ws2, src );
    const PhiVolT s3( ws3, src );
    
    SpatFldPtr<PhiVolT> fdir = SpatialFieldStore::get<PhiVolT>( aVel );
    *fdir <<= - pm[direc] * aVel; // flow direction

    SpatFldPtr<PhiVolT> r = SpatialFieldStore::get<PhiVolT>( s1 );
    *r <<= (s3-s2)/(s2-s1);
    
    const bool isBoundaryFace = (hasMinusBoundary_ && direc==0) || (hasPlusBoundary_ && direc==1);    

    switch (limiterType_) {
    case Wasatch::UPWIND  : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, UPWIND   ); break;
    case Wasatch::SUPERBEE: d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, SUPERBEE ); break;
    case Wasatch::CHARM   : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, CHARM    ); break;
    case Wasatch::KOREN   : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, KOREN    ); break;
    case Wasatch::MC      : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, MC       ); break;
    case Wasatch::OSPRE   : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, OSPRE    ); break;
    case Wasatch::SMART   : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, SMART    ); break;
    case Wasatch::VANLEER : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, VANLEER  ); break;
    case Wasatch::HCUS    : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, HCUS     ); break;
    case Wasatch::MINMOD  : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, MINMOD   ); break;
    case Wasatch::HQUICK  : d <<= CALCULATE_BOUNDARY_LIMITER(*fdir, *r, HQUICK   ); break;
    default               : d <<= 0.0;                                              break;
    }            
  }

  // now do interior
  destExtent = wdest.extent() - unitNormal_*3 - wdest.has_bc()*unitNormal_;
  destBaseOffset = wdest.offset() + unitNormal_*2;
  const MemoryWindow wd( wdest.glob_dim(),
                         destBaseOffset,
                         destExtent,
                         wdest.has_bc(0), wdest.has_bc(1), wdest.has_bc(2) );

  PhiVolT d( wd, destVals, ExternalStorage, dMemType, dDevIdx );
  const PhiVolT aVel( wd, velVals, ExternalStorage, advelMemType, advelDevIdx );
  
  // build the source fields - these correspond to windows into minus-minus,
  // minus, plus, and plus-plus with respect to destination (face).
  build_src_fields(src);
  
  const PhiVolT& smm = *srcFields_[0];
  const PhiVolT& sm  = *srcFields_[1];
  const PhiVolT& sp  = *srcFields_[2];
  const PhiVolT& spp = *srcFields_[3];
  
  SpatFldPtr<PhiVolT> rm = SpatialFieldStore::get<PhiVolT>( sm );
  *rm <<= (sm-smm)/(sp-sm);
  
  SpatFldPtr<PhiVolT> rp = SpatialFieldStore::get<PhiVolT>( sp );
  *rp <<= (spp-sp)/(sp-sm);
      
  switch (limiterType_) {
  case Wasatch::UPWIND  : d <<= 0.0;                                                   break;
  case Wasatch::SUPERBEE: d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, SUPERBEE ); break;
  case Wasatch::CHARM   : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, CHARM    ); break;
  case Wasatch::KOREN   : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, KOREN    ); break;
  case Wasatch::MC      : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, MC       ); break;
  case Wasatch::OSPRE   : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, OSPRE    ); break;
  case Wasatch::SMART   : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, SMART    ); break;
  case Wasatch::VANLEER : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, VANLEER  ); break;
  case Wasatch::HCUS    : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, HCUS     ); break;
  case Wasatch::MINMOD  : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, MINMOD   ); break;
  case Wasatch::HQUICK  : d <<= CALCULATE_INTERIOR_LIMITER(aVel, *rm, *rp, HQUICK   ); break;
  default               : d <<= 0.0;                                                   break;
  }

# ifdef ENABLE_THREADS
  mutex_.unlock();
# endif
}

//==================================================================
// Explicit template instantiation
namespace SS = SpatialOps::structured;

template class FluxLimiterInterpolant< SS::SVolField, SS::SSurfXField >;
template class FluxLimiterInterpolant< SS::SVolField, SS::SSurfYField >;
template class FluxLimiterInterpolant< SS::SVolField, SS::SSurfZField >;

//==================================================================
