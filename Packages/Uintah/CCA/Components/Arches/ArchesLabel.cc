//----- ArchesLabel.cc ----------------------------------------------

#include <Packages/Uintah/CCA/Components/Arches/ArchesLabel.h>
#include <Packages/Uintah/CCA/Components/Arches/CellInformationP.h>
#include <Packages/Uintah/CCA/Components/Arches/CellInformation.h>
#include <Packages/Uintah/Core/Grid/VarLabel.h>
#include <Packages/Uintah/Core/Grid/PerPatch.h>
#include <Packages/Uintah/Core/Grid/CCVariable.h>
#include <Packages/Uintah/Core/Grid/SFCXVariable.h>
#include <Packages/Uintah/Core/Grid/SFCYVariable.h>
#include <Packages/Uintah/Core/Grid/SFCZVariable.h>
#include <Packages/Uintah/Core/Grid/VarTypes.h>
#include <Packages/Uintah/Core/Grid/SimulationState.h>

using namespace Uintah;

//****************************************************************************
// Default constructor for ArchesLabel
//****************************************************************************
ArchesLabel::ArchesLabel()
{
  int noStencil = 7;
  d_stencilMatl = scinew MaterialSubset();
  for (int i = 0; i < noStencil; i++)
    d_stencilMatl->add(i);
  d_stencilMatl->addReference();

  // Cell Information
  d_cellInfoLabel = VarLabel::create("cellInformation",
			    PerPatch<CellInformationP>::getTypeDescription());
  // Cell type
  d_cellTypeLabel = VarLabel::create("cellType", 
				  CCVariable<int>::getTypeDescription() );
  // labels for inlet and outlet flow rates
  d_totalflowINLabel = VarLabel::create("totalflowIN",
     ReductionVariable<double, Reductions::Sum<double> >::getTypeDescription()); 
  d_totalflowOUTLabel = VarLabel::create("totalflowOUT",
     ReductionVariable<double, Reductions::Sum<double> >::getTypeDescription()); 
  d_totalflowOUToutbcLabel = VarLabel::create("totalflowOUToutbc",
     ReductionVariable<double, Reductions::Sum<double> >::getTypeDescription()); 
  d_totalAreaOUTLabel = VarLabel::create("totalAreaOUT",
     ReductionVariable<double, Reductions::Sum<double> >::getTypeDescription()); 
  d_denAccumLabel = VarLabel::create("denAccum",
     ReductionVariable<double, Reductions::Sum<double> >::getTypeDescription()); 

  // Density Labels
  d_densityINLabel = VarLabel::create("densityIN", 
				   CCVariable<double>::getTypeDescription() );
  d_densitySPLabel = VarLabel::create("densitySP", 
				   CCVariable<double>::getTypeDescription() );
  d_densityCPLabel = VarLabel::create("densityCP", 
				  CCVariable<double>::getTypeDescription() );
  // Viscosity Labels
  d_viscosityINLabel = VarLabel::create("viscosityIN", 
				   CCVariable<double>::getTypeDescription() );
  d_viscosityCTSLabel = VarLabel::create("viscosityCTS", 
				      CCVariable<double>::getTypeDescription() );
  // Pressure Labels
  d_pressureINLabel = VarLabel::create("pressureIN", 
				    CCVariable<double>::getTypeDescription() );
  d_pressureSPBCLabel = VarLabel::create("pressureSPBC", 
				      CCVariable<double>::getTypeDescription() );
  d_pressurePSLabel = VarLabel::create("pressurePS", 
				    CCVariable<double>::getTypeDescription() );
  // Pressure Coeff Labels
  d_presCoefPBLMLabel = VarLabel::create("presCoefPBLM", 
				      CCVariable<double>::getTypeDescription() );
  // Pressure Linear Src Labels
  d_presLinSrcPBLMLabel = VarLabel::create("presLinSrcPBLM", 
					CCVariable<double>::getTypeDescription() );
  // Pressure Non Linear Src Labels
  d_presNonLinSrcPBLMLabel = VarLabel::create("presNonLinSrcPBLM", 
					   CCVariable<double>::getTypeDescription() );
  // U-Velocity Labels
  d_uVelocityINLabel = VarLabel::create("uVelocityIN", 
				    SFCXVariable<double>::getTypeDescription() );
  d_uVelocitySPLabel = VarLabel::create("uVelocitySP", 
				    SFCXVariable<double>::getTypeDescription() );
  d_uVelocitySPBCLabel = VarLabel::create("uVelocitySPBC", 
				       SFCXVariable<double>::getTypeDescription() );
  d_uVelocitySIVBCLabel = VarLabel::create("uVelocitySIVBC", 
				       SFCXVariable<double>::getTypeDescription() );
  d_uVelocityCPBCLabel = VarLabel::create("uVelocityCPBC", 
				       SFCXVariable<double>::getTypeDescription() );
  // U-Velocity Coeff Labels
  d_uVelCoefPBLMLabel = VarLabel::create("uVelCoefPBLM",
			       SFCXVariable<double>::getTypeDescription() );
  // U-Velocity Convection Coeff Labels
  d_uVelConvCoefPBLMLabel = VarLabel::create("uVelConvCoefPBLM",
				        SFCXVariable<double>::getTypeDescription() );
  // U-Velocity Linear Src Labels
  d_uVelLinSrcPBLMLabel = VarLabel::create("uVelLinSrcPBLM",
				 SFCXVariable<double>::getTypeDescription() );
  // U-Velocity Non Linear Src Labels
  d_uVelNonLinSrcPBLMLabel = VarLabel::create("uVelNonLinSrcPBLM",
				    SFCXVariable<double>::getTypeDescription() );
  // U-Velocity Coeff Labels
  d_uVelCoefMBLMLabel = VarLabel::create("uVelCoefMBLM",
			       SFCXVariable<double>::getTypeDescription() );
  // U-Velocity Convection Coeff Labels
  d_uVelConvCoefMBLMLabel = VarLabel::create("uVelConvCoefMBLM",
				        SFCXVariable<double>::getTypeDescription() );
  // U-Velocity Linear Src Labels
  d_uVelLinSrcMBLMLabel = VarLabel::create("uVelLinSrcMBLM",
				 SFCXVariable<double>::getTypeDescription() );
  // U-Velocity Non Linear Src Labels
  d_uVelNonLinSrcMBLMLabel = VarLabel::create("uVelNonLinSrcMBLM",
				    SFCXVariable<double>::getTypeDescription() );
  // V-Velocity Labels
  d_vVelocityINLabel = VarLabel::create("vVelocityIN", 
				    SFCYVariable<double>::getTypeDescription() );
  d_vVelocitySPLabel = VarLabel::create("vVelocitySP", 
				    SFCYVariable<double>::getTypeDescription() );
  d_vVelocitySPBCLabel = VarLabel::create("vVelocitySPBC", 
				       SFCYVariable<double>::getTypeDescription() );
  d_vVelocitySIVBCLabel = VarLabel::create("vVelocitySIVBC", 
				       SFCYVariable<double>::getTypeDescription() );
  d_vVelocityCPBCLabel = VarLabel::create("vVelocityCPBC", 
				       SFCYVariable<double>::getTypeDescription() );
  // V-Velocity Coeff Labels
  d_vVelCoefPBLMLabel = VarLabel::create("vVelCoefPBLM",
			       SFCYVariable<double>::getTypeDescription() );
  // V-Velocity Convection Coeff Labels
  d_vVelConvCoefPBLMLabel = VarLabel::create("vVelConvCoefPBLM",
				   SFCYVariable<double>::getTypeDescription() );
  // V-Velocity Linear Src Labels
  d_vVelLinSrcPBLMLabel = VarLabel::create("vVelLinSrcPBLM",
				 SFCYVariable<double>::getTypeDescription() );
  // V-Velocity Non Linear Src Labels
  d_vVelNonLinSrcPBLMLabel = VarLabel::create("vVelNonLinSrcPBLM",
				    SFCYVariable<double>::getTypeDescription() );
  // V-Velocity Coeff Labels
  d_vVelCoefMBLMLabel = VarLabel::create("vVelCoefMBLM",
			       SFCYVariable<double>::getTypeDescription() );
  // V-Velocity Convection Coeff Labels
  d_vVelConvCoefMBLMLabel = VarLabel::create("vVelConvCoefMBLM",
				   SFCYVariable<double>::getTypeDescription() );
  // V-Velocity Linear Src Labels
  d_vVelLinSrcMBLMLabel = VarLabel::create("vVelLinSrcMBLM",
				 SFCYVariable<double>::getTypeDescription() );

  // labels for ref density and pressure
  d_refDensity_label = VarLabel::create("refDensityLabel",
				       sum_vartype::getTypeDescription() );
  d_refPressure_label = VarLabel::create("refPressureLabel",
				       sum_vartype::getTypeDescription() );

  // V-Velocity Non Linear Src Labels
  d_vVelNonLinSrcMBLMLabel = VarLabel::create("vVelNonLinSrcMBLM",
				    SFCYVariable<double>::getTypeDescription() );
  // W-Velocity Labels
  d_wVelocityINLabel = VarLabel::create("wVelocityIN", 
				    SFCZVariable<double>::getTypeDescription() );
  d_wVelocitySPLabel = VarLabel::create("wVelocitySP", 
				    SFCZVariable<double>::getTypeDescription() );
  d_wVelocitySPBCLabel = VarLabel::create("wVelocitySPBC", 
				       SFCZVariable<double>::getTypeDescription() );
  d_wVelocitySIVBCLabel = VarLabel::create("wVelocitySIVBC", 
				       SFCZVariable<double>::getTypeDescription() );
  d_wVelocityCPBCLabel = VarLabel::create("wVelocityCPBC", 
				       SFCZVariable<double>::getTypeDescription() );
  // W-Velocity Coeff Labels
  d_wVelCoefPBLMLabel = VarLabel::create("wVelCoefPBLM",
			       SFCZVariable<double>::getTypeDescription() );
  // W-Velocity Convection Coeff Labels
  d_wVelConvCoefPBLMLabel = VarLabel::create("wVelConvCoefPBLM",
				   SFCZVariable<double>::getTypeDescription() );
  // W-Velocity Linear Src Labels
  d_wVelLinSrcPBLMLabel = VarLabel::create("wVelLinSrcPBLM",
				 SFCZVariable<double>::getTypeDescription() );
  // W-Velocity Non Linear Src Labels
  d_wVelNonLinSrcPBLMLabel = VarLabel::create("wVelNonLinSrcPBLM",
				    SFCZVariable<double>::getTypeDescription() );
  // W-Velocity Coeff Labels
  d_wVelCoefMBLMLabel = VarLabel::create("wVelCoefMBLM",
			       SFCZVariable<double>::getTypeDescription() );
  // W-Velocity Convection Coeff Labels
  d_wVelConvCoefMBLMLabel = VarLabel::create("wVelConvCoefMBLM",
				   SFCZVariable<double>::getTypeDescription() );
  // W-Velocity Linear Src Labels
  d_wVelLinSrcMBLMLabel = VarLabel::create("wVelLinSrcMBLM",
				 SFCZVariable<double>::getTypeDescription() );
  // W-Velocity Non Linear Src Labels
  d_wVelNonLinSrcMBLMLabel = VarLabel::create("wVelNonLinSrcMBLM",
				    SFCZVariable<double>::getTypeDescription() );
  // Scalar 
  d_scalarINLabel = VarLabel::create("scalarIN",
				    CCVariable<double>::getTypeDescription() );
  d_scalarSPLabel = VarLabel::create("scalarSP",
				   CCVariable<double>::getTypeDescription() );

  d_scalarCPBCLabel = VarLabel::create("scalarCPBC",
				      CCVariable<double>::getTypeDescription() );

  // scalar variance
  d_scalarVarINLabel = VarLabel::create("scalarVarIN", 
				       CCVariable<double>::getTypeDescription() );

  d_scalarVarSPLabel = VarLabel::create("scalarVarSP", 
				       CCVariable<double>::getTypeDescription() );
  // Scalar Coef
  d_scalCoefSBLMLabel = VarLabel::create("scalCoefSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Conv Coef
  d_scalConvCoefSBLMLabel = VarLabel::create("scalConvCoefSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Linear Src
  d_scalLinSrcSBLMLabel = VarLabel::create("scalLinSrcSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Non Linear Src
  d_scalNonLinSrcSBLMLabel = VarLabel::create("scalNonLinSrcSBLM",
				   CCVariable<double>::getTypeDescription() );

  // reactive scalar

  d_reactscalarINLabel = VarLabel::create("reactscalarIN",
				    CCVariable<double>::getTypeDescription() );
  d_reactscalarSPLabel = VarLabel::create("reactscalarSP",
				   CCVariable<double>::getTypeDescription() );

  d_reactscalarOUTBCLabel = VarLabel::create("reactscalarOUTBC",
				   CCVariable<double>::getTypeDescription() );

  d_reactscalarCPBCLabel = VarLabel::create("reactscalarCPBC",
				      CCVariable<double>::getTypeDescription() );

  // reactscalar variance
  d_reactscalarVarINLabel = VarLabel::create("reactscalarVarIN", 
				       CCVariable<double>::getTypeDescription() );

  d_reactscalarVarSPLabel = VarLabel::create("reactscalarVarSP", 
				       CCVariable<double>::getTypeDescription() );
  // Reactscalar Coef
  d_reactscalCoefSBLMLabel = VarLabel::create("reactscalCoefSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Conv Coef
  d_reactscalConvCoefSBLMLabel = VarLabel::create("reactscalConvCoefSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Linear Src
  d_reactscalLinSrcSBLMLabel = VarLabel::create("reactscalLinSrcSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Non Linear Src
  d_reactscalNonLinSrcSBLMLabel = VarLabel::create("reactscalNonLinSrcSBLM",
				   CCVariable<double>::getTypeDescription() );


  // labels for nonlinear residuals
  d_presResidPSLabel = VarLabel::create("presResidPSLabel",
				        ReductionVariable<double,
				       Reductions::Sum<double> >::getTypeDescription());
  d_presTruncPSLabel = VarLabel::create("presTruncPSLabel",
				       sum_vartype::getTypeDescription() );
  d_uVelResidPSLabel = VarLabel::create("uVelResidPSLabel",
				       sum_vartype::getTypeDescription() );
  d_uVelTruncPSLabel = VarLabel::create("uVelTruncPSLabel",
				       sum_vartype::getTypeDescription() );
  d_vVelResidPSLabel = VarLabel::create("vVelResidPSLabel",
				       sum_vartype::getTypeDescription() );
  d_vVelTruncPSLabel = VarLabel::create("vVelTruncPSLabel",
				       sum_vartype::getTypeDescription() );
  d_wVelResidPSLabel = VarLabel::create("wVelResidPSLabel",
				       sum_vartype::getTypeDescription() );
  d_wVelTruncPSLabel = VarLabel::create("wVelTruncPSLabel",
				       sum_vartype::getTypeDescription() );
  d_scalarResidLabel = VarLabel::create("scalarResidLabel",
				       sum_vartype::getTypeDescription() );
  d_scalarTruncLabel = VarLabel::create("scalarTruncLabel",
				       sum_vartype::getTypeDescription() );


  d_pressureRes = VarLabel::create("pressureRes",
				   CCVariable<double>::getTypeDescription() );
  d_uVelocityRes = VarLabel::create("uVelocityRes",
				   SFCXVariable<double>::getTypeDescription() );
  d_vVelocityRes = VarLabel::create("vVelocityRes",
				   SFCYVariable<double>::getTypeDescription() );
  d_wVelocityRes = VarLabel::create("wVelocityRes",
				   SFCZVariable<double>::getTypeDescription() );
  d_scalarRes = VarLabel::create("scalarRes",
				   CCVariable<double>::getTypeDescription() );

  d_reactscalarRes = VarLabel::create("reactscalarRes",
				   CCVariable<double>::getTypeDescription() );

  d_enthalpyRes = VarLabel::create("enthalpyRes",
				   CCVariable<double>::getTypeDescription() );

  d_old_uVelocityGuess = VarLabel::create("olduVelocityguess",
				       SFCXVariable<double>::getTypeDescription());
  d_old_vVelocityGuess = VarLabel::create("oldvVelocityguess",
				       SFCYVariable<double>::getTypeDescription());
  d_old_wVelocityGuess = VarLabel::create("oldwVelocityguess",
				       SFCZVariable<double>::getTypeDescription());
  d_old_scalarGuess = VarLabel::create("oldscalarguess",
				       CCVariable<double>::getTypeDescription());

  // Unsure stuff
  // Unsure stuff
  d_DUPBLMLabel = VarLabel::create("DUPBLM",
				SFCXVariable<double>::getTypeDescription() );
  d_DVPBLMLabel = VarLabel::create("DVPBLM",
				SFCYVariable<double>::getTypeDescription() );
  d_DWPBLMLabel = VarLabel::create("DWPBLM",
				SFCZVariable<double>::getTypeDescription() );
  d_DUMBLMLabel = VarLabel::create("DUMBLM",
				SFCXVariable<double>::getTypeDescription() );
  d_DVMBLMLabel = VarLabel::create("DVMBLM",
				SFCYVariable<double>::getTypeDescription() );
  d_DWMBLMLabel = VarLabel::create("DWMBLM",
				SFCZVariable<double>::getTypeDescription() );

  // Labels that access the velocity stored as a cell centered vector
  // after interpolation (for use in visualization)

  d_oldCCVelocityLabel = VarLabel::create("oldCCVelocity",
				CCVariable<Vector>::getTypeDescription() );
  d_newCCVelocityLabel = VarLabel::create("newCCVelocity",
				CCVariable<Vector>::getTypeDescription() );
  d_newCCUVelocityLabel = VarLabel::create("newCCUVelocity",
				CCVariable<double>::getTypeDescription() );
  d_newCCVVelocityLabel = VarLabel::create("newCCVVelocity",
				CCVariable<double>::getTypeDescription() );
  d_newCCWVelocityLabel = VarLabel::create("newCCWVelocity",
				CCVariable<double>::getTypeDescription() );

  // for pressure grad term in momentum

  d_pressGradUSuLabel = VarLabel::create("pressGradUSu",
					SFCXVariable<double>::getTypeDescription() );
  d_pressGradVSuLabel = VarLabel::create("pressGradVSu",
					SFCYVariable<double>::getTypeDescription() );
  d_pressGradWSuLabel = VarLabel::create("pressGradWSu",
					SFCZVariable<double>::getTypeDescription() );
  // multimaterial labels

  // multimaterial wall/intrusion cells

  d_mmcellTypeLabel = VarLabel::create("mmcellType",
				      CCVariable<int>::getTypeDescription() );

  // Label for void fraction, after correction for wall cells using cutoff

  d_mmgasVolFracLabel = VarLabel::create("mmgasVolFrac",
					CCVariable<double>::getTypeDescription() );

    // for reacting flows
  d_tempINLabel = VarLabel::create("tempIN",
				  CCVariable<double>::getTypeDescription() );
  d_co2INLabel = VarLabel::create("co2IN",
				  CCVariable<double>::getTypeDescription() );


  // Array containing the reference density multiplied by the void fraction
  // used for correct reference density subtraction in the multimaterial
  // case

  d_denRefArrayLabel = VarLabel::create("denRefArray",
					CCVariable<double>::getTypeDescription() );

  // Microscopic density (i.e., without void fraction) of gas

  d_densityMicroLabel = VarLabel::create("denMicro",
					CCVariable<double>::getTypeDescription() );
  d_densityMicroINLabel = VarLabel::create("denMicroIN",
					CCVariable<double>::getTypeDescription() );

  // Label for the sum of the relative pressure and the hydrostatic 
  // contribution

  d_pressPlusHydroLabel = VarLabel::create("pPlusHydro",
					CCVariable<double>::getTypeDescription() );

  d_scalarOUTBCLabel =  VarLabel::create("scalarOUTBC",
					CCVariable<double>::getTypeDescription() );
  d_uVelocityOUTBCLabel = VarLabel::create("uVelocityOUTBC",
				   SFCXVariable<double>::getTypeDescription() );
  d_vVelocityOUTBCLabel= VarLabel::create("vVelocityOUTBC",
				  SFCYVariable<double>::getTypeDescription() );
  d_wVelocityOUTBCLabel= VarLabel::create("wVelocityOUTBC",
				   SFCZVariable<double>::getTypeDescription() );

  d_uvwoutLabel = VarLabel::create("uvwout",
	  ReductionVariable<double, Reductions::Min<double> >::getTypeDescription()); 


  // predictor-corrector labels
  // Scalar Coef
  d_scalCoefPredLabel = VarLabel::create("scalCoefPred",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Conv Coef
  d_scalConvCoefPredLabel = VarLabel::create("scalConvCoefPred",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Linear Src
  d_scalLinSrcPredLabel = VarLabel::create("scalLinSrcPred",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Non Linear Src
  d_scalNonLinSrcPredLabel = VarLabel::create("scalNonLinSrcPred",
				   CCVariable<double>::getTypeDescription() );

  d_scalCoefCorrLabel = VarLabel::create("scalCoefCorr",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Conv Coef
  d_scalConvCoefCorrLabel = VarLabel::create("scalConvCoefCorr",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Linear Src
  d_scalLinSrcCorrLabel = VarLabel::create("scalLinSrcCorr",
				   CCVariable<double>::getTypeDescription() );
  // Scalar Non Linear Src
  d_scalNonLinSrcCorrLabel = VarLabel::create("scalNonLinSrcCorr",
				   CCVariable<double>::getTypeDescription() );

  d_scalarPredLabel = VarLabel::create("scalarPred",
				   CCVariable<double>::getTypeDescription() );


  // predictor-corrector labels
  // Reactscalar Coef
  d_reactscalCoefPredLabel = VarLabel::create("reactscalCoefPred",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Conv Coef
  d_reactscalConvCoefPredLabel = VarLabel::create("reactscalConvCoefPred",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Linear Src
  d_reactscalLinSrcPredLabel = VarLabel::create("reactscalLinSrcPred",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Non Linear Src
  d_reactscalNonLinSrcPredLabel = VarLabel::create("reactscalNonLinSrcPred",
				   CCVariable<double>::getTypeDescription() );

  d_reactscalCoefCorrLabel = VarLabel::create("reactscalCoefCorr",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Conv Coef
  d_reactscalConvCoefCorrLabel = VarLabel::create("reactscalConvCoefCorr",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Linear Src
  d_reactscalLinSrcCorrLabel = VarLabel::create("reactscalLinSrcCorr",
				   CCVariable<double>::getTypeDescription() );
  // Reactscalar Non Linear Src
  d_reactscalNonLinSrcCorrLabel = VarLabel::create("reactscalNonLinSrcCorr",
				   CCVariable<double>::getTypeDescription() );

  d_reactscalarPredLabel = VarLabel::create("reactscalarPred",
				   CCVariable<double>::getTypeDescription() );

  d_densityPredLabel = VarLabel::create("densityPred",
				   CCVariable<double>::getTypeDescription() );

  d_uVelRhoHatLabel = VarLabel::create("uvelRhoHat",
				   SFCXVariable<double>::getTypeDescription() );
  d_vVelRhoHatLabel= VarLabel::create("vvelRhoHat",
				  SFCYVariable<double>::getTypeDescription() );
  d_wVelRhoHatLabel= VarLabel::create("wvelRhoHat",
				   SFCZVariable<double>::getTypeDescription() );

  d_pressurePredLabel = VarLabel::create("pressurePred", 
				    CCVariable<double>::getTypeDescription() );
  // Pressure Coeff Labels
  d_presCoefCorrLabel = VarLabel::create("presCoefCorr", 
				      CCVariable<double>::getTypeDescription() );
  // Pressure Linear Src Labels
  d_presLinSrcCorrLabel = VarLabel::create("presLinSrcCorr", 
					CCVariable<double>::getTypeDescription() );
  // Pressure Non Linear Src Labels
  d_presNonLinSrcCorrLabel = VarLabel::create("presNonLinSrcCorr", 
					   CCVariable<double>::getTypeDescription() );
  // U-Velocity Labels
  d_uVelocityPredLabel = VarLabel::create("uVelocityPred", 
				    SFCXVariable<double>::getTypeDescription() );
  d_vVelocityPredLabel = VarLabel::create("vVelocityPred", 
				    SFCYVariable<double>::getTypeDescription() );
  d_wVelocityPredLabel = VarLabel::create("wVelocityPred", 
				    SFCZVariable<double>::getTypeDescription() );

  // enthalpy labels
  // Enthalpy 
  d_enthalpyINLabel = VarLabel::create("enthalpyIN",
				    CCVariable<double>::getTypeDescription() );
  d_enthalpySPLabel = VarLabel::create("enthalpySP",
				   CCVariable<double>::getTypeDescription() );
  d_enthalpySPBCLabel = VarLabel::create("enthalpySPBC",
				   CCVariable<double>::getTypeDescription() );

  d_enthalpyRXNLabel = VarLabel::create("enthalpyRXN",
				   CCVariable<double>::getTypeDescription() );

  d_enthalpyCPBCLabel = VarLabel::create("enthalpyCPBC",
				      CCVariable<double>::getTypeDescription() );
  d_enthalpyOUTBCLabel = VarLabel::create("enthalpyOUTBC",
				      CCVariable<double>::getTypeDescription() );

  // Enthalpy Coef
  d_enthCoefSBLMLabel = VarLabel::create("enthCoefSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Conv Coef
  d_enthConvCoefSBLMLabel = VarLabel::create("enthConvCoefSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Linear Src
  d_enthLinSrcSBLMLabel = VarLabel::create("enthLinSrcSBLM",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Non Linear Src
  d_enthNonLinSrcSBLMLabel = VarLabel::create("enthNonLinSrcSBLM",
				   CCVariable<double>::getTypeDescription() );

  // predictor-corrector labels
  // Enthalpy Coef
  d_enthCoefPredLabel = VarLabel::create("enthCoefPred",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Conv Coef
  d_enthConvCoefPredLabel = VarLabel::create("enthConvCoefPred",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Linear Src
  d_enthLinSrcPredLabel = VarLabel::create("enthLinSrcPred",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Non Linear Src
  d_enthNonLinSrcPredLabel = VarLabel::create("enthNonLinSrcPred",
				   CCVariable<double>::getTypeDescription() );

  d_enthCoefCorrLabel = VarLabel::create("enthCoefCorr",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Conv Coef
  d_enthConvCoefCorrLabel = VarLabel::create("enthConvCoefCorr",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Linear Src
  d_enthLinSrcCorrLabel = VarLabel::create("enthLinSrcCorr",
				   CCVariable<double>::getTypeDescription() );
  // Enthalpy Non Linear Src
  d_enthNonLinSrcCorrLabel = VarLabel::create("enthNonLinSrcCorr",
				   CCVariable<double>::getTypeDescription() );

  d_enthalpyPredLabel = VarLabel::create("enthalpyPred",
				   CCVariable<double>::getTypeDescription() );
  // for radiation
  d_absorpINLabel = VarLabel::create("absorpIN",
				    CCVariable<double>::getTypeDescription() );
  d_sootFVINLabel = VarLabel::create("sootFVIN",
				    CCVariable<double>::getTypeDescription() );

  d_reactscalarSRCINLabel = VarLabel::create("reactscalarSRCIN",
				    CCVariable<double>::getTypeDescription() );


}

//****************************************************************************
// Destructor
//****************************************************************************
ArchesLabel::~ArchesLabel()
{
  if (d_stencilMatl->removeReference())
    delete d_stencilMatl;
}

void ArchesLabel::setSharedState(SimulationStateP& sharedState)
{
  d_sharedState = sharedState;
}
