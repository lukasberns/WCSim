// modified from Geant4's processes/optical/src/G4OpBoundaryProcess.cc 
//
// Author:      Lukas Berns
// mail:        lukas.berns@hep.phys.titech.ac.jp
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpProcessSubType.hh"

#include "WCSimOpRetroBoundaryProcess.hh"
#include "G4GeometryTolerance.hh"

#include "G4VSensitiveDetector.hh"
#include "G4ParallelWorldProcess.hh"

#include "G4SystemOfUnits.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

// WCSimOpRetroBoundaryProcess::operator=(const G4OpBoundaryProcess &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

WCSimOpRetroBoundaryProcess::WCSimOpRetroBoundaryProcess(const G4String& processName,
                                               G4ProcessType type)
             : G4VDiscreteProcess(processName, type)
{
        if ( verboseLevel > 0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

        SetProcessSubType(fOpBoundary);

	theStatus = Undefined;
	theModel = prism;
	theFinish = perfect;
        theReflectivity =  1.;
        theEfficiency   =  0.;
        theTransmittance = 0.;

        theSurfaceRoughness = 0.;

        prob_sl = 0.;
        prob_ss = 0.;
        prob_bs = 0.;

        PropertyPointer  = NULL;
        PropertyPointer1 = NULL;
        PropertyPointer2 = NULL;

        Material1 = NULL;
        Material2 = NULL;

        OpticalSurface = NULL;

        kCarTolerance = G4GeometryTolerance::GetInstance()
                        ->GetSurfaceTolerance();

        iTE = iTM = 0;
        thePhotonMomentum = 0.;
        Rindex1 = Rindex2 = cost1 = cost2 = sint1 = sint2 = 0.;

        idx = idy = 0;
        DichroicVector = NULL;
}

// WCSimOpRetroBoundaryProcess::WCSimOpRetroBoundaryProcess(const G4OpBoundaryProcess &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

WCSimOpRetroBoundaryProcess::~WCSimOpRetroBoundaryProcess(){}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// ------------
//

G4VParticleChange*
WCSimOpRetroBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
        theStatus = Undefined;

        aParticleChange.Initialize(aTrack);
        aParticleChange.ProposeVelocity(aTrack.GetVelocity());

        // Get hyperStep from  G4ParallelWorldProcess
        //  NOTE: PostSetpDoIt of this process should be
        //        invoked after G4ParallelWorldProcess!

        const G4Step* pStep = &aStep;

        const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();
        
        if (hStep) pStep = hStep;

        G4bool isOnBoundary =
                (pStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);

        if (isOnBoundary) {
           Material1 = pStep->GetPreStepPoint()->GetMaterial();
           Material2 = pStep->GetPostStepPoint()->GetMaterial();
        } else {
           theStatus = NotAtBoundary;
           if ( verboseLevel > 0) BoundaryProcessVerbose();
           return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
        }

        G4VPhysicalVolume* thePrePV  =
                               pStep->GetPreStepPoint() ->GetPhysicalVolume();
        G4VPhysicalVolume* thePostPV =
                               pStep->GetPostStepPoint()->GetPhysicalVolume();

        if ( verboseLevel > 0 ) {
           G4cout << " Photon at Boundary! " << G4endl;
           if (thePrePV)  G4cout << " thePrePV:  " << thePrePV->GetName()  << G4endl;
           if (thePostPV) G4cout << " thePostPV: " << thePostPV->GetName() << G4endl;
        }

	if (aTrack.GetStepLength()<=kCarTolerance/2){
	        theStatus = StepTooSmall;
                if ( verboseLevel > 0) BoundaryProcessVerbose();
	        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

	thePhotonMomentum = aParticle->GetTotalMomentum();
        OldMomentum       = aParticle->GetMomentumDirection();
	OldPolarization   = aParticle->GetPolarization();

        if ( verboseLevel > 0 ) {
           G4cout << " Old Momentum Direction: " << OldMomentum     << G4endl;
           G4cout << " Old Polarization:       " << OldPolarization << G4endl;
        }

        G4ThreeVector theGlobalPoint = pStep->GetPostStepPoint()->GetPosition();

        G4bool valid;
        //  Use the new method for Exit Normal in global coordinates,
        //    which provides the normal more reliably.

        // ID of Navigator which limits step

        G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
        std::vector<G4Navigator*>::iterator iNav =
                G4TransportationManager::GetTransportationManager()->
                                         GetActiveNavigatorsIterator();
        
        theGlobalNormal =
                   (iNav[hNavId])->GetGlobalExitNormal(theGlobalPoint,&valid);

        if (valid) {
          theGlobalNormal = -theGlobalNormal;
        }
        else 
        {
          G4ExceptionDescription ed;
          ed << " WCSimOpRetroBoundaryProcess/PostStepDoIt(): "
                 << " The Navigator reports that it returned an invalid normal"
                 << G4endl;
          G4Exception("WCSimOpRetroBoundaryProcess::PostStepDoIt", "OpBoun01",
                      EventMustBeAborted,ed,
                      "Invalid Surface Normal - Geometry must return valid surface normal");
        }

        if (OldMomentum * theGlobalNormal > 0.0) {
#ifdef G4OPTICAL_DEBUG
           G4ExceptionDescription ed;
           ed << " WCSimOpRetroBoundaryProcess/PostStepDoIt(): "
              << " theGlobalNormal points in a wrong direction. "
              << G4endl;
           ed << "    The momentum of the photon arriving at interface (oldMomentum)"
              << " must exit the volume cross in the step. " << G4endl;
           ed << "  So it MUST have dot < 0 with the normal that Exits the new volume (globalNormal)." << G4endl;
           ed << "  >> The dot product of oldMomentum and global Normal is " << OldMomentum*theGlobalNormal << G4endl;
           ed << "     Old Momentum  (during step)     = " << OldMomentum << G4endl;
           ed << "     Global Normal (Exiting New Vol) = " << theGlobalNormal << G4endl;
           ed << G4endl;
           G4Exception("WCSimOpRetroBoundaryProcess::PostStepDoIt", "OpBoun02",
                       EventMustBeAborted,  // Or JustWarning to see if it happens repeatedbly on one ray
                       ed,
                      "Invalid Surface Normal - Geometry must return valid surface normal pointing in the right direction");
#else
           theGlobalNormal = -theGlobalNormal;
#endif
        }

	G4MaterialPropertiesTable* aMaterialPropertiesTable;
        G4MaterialPropertyVector* Rindex;

	aMaterialPropertiesTable = Material1->GetMaterialPropertiesTable();
        if (aMaterialPropertiesTable) {
		Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	}
	else {
	        theStatus = NoRINDEX;
                if ( verboseLevel > 0) BoundaryProcessVerbose();
                aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        if (Rindex) {
		Rindex1 = Rindex->Value(thePhotonMomentum);
	}
	else {
	        theStatus = NoRINDEX;
                if ( verboseLevel > 0) BoundaryProcessVerbose();
                aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        theReflectivity =  1.;
        theEfficiency   =  0.;
        theTransmittance = 0.;

        theSurfaceRoughness = 0.;

        theModel = prism;
        theFinish = perfect;

        G4SurfaceType type = dielectric_dielectric;

        Rindex = NULL;
        OpticalSurface = NULL;

        G4LogicalSurface* Surface = NULL;

        Surface = G4LogicalBorderSurface::GetSurface(thePrePV, thePostPV);

        if (Surface == NULL){
          G4bool enteredDaughter= (thePostPV->GetMotherLogical() ==
                                   thePrePV ->GetLogicalVolume());
	  if(enteredDaughter){
	    Surface = 
              G4LogicalSkinSurface::GetSurface(thePostPV->GetLogicalVolume());
	    if(Surface == NULL)
	      Surface =
                G4LogicalSkinSurface::GetSurface(thePrePV->GetLogicalVolume());
	  }
	  else {
	    Surface =
              G4LogicalSkinSurface::GetSurface(thePrePV->GetLogicalVolume());
	    if(Surface == NULL)
	      Surface =
                G4LogicalSkinSurface::GetSurface(thePostPV->GetLogicalVolume());
	  }
	}

	if (Surface) OpticalSurface = 
           dynamic_cast <WCSimOpticalSurface*> (Surface->GetSurfaceProperty());

	if (OpticalSurface) {

           type      = OpticalSurface->GetType();
	   theModel  = OpticalSurface->GetModel();
	   theFinish = OpticalSurface->GetFinish();

	   aMaterialPropertiesTable = OpticalSurface->
					GetMaterialPropertiesTable();

           if (aMaterialPropertiesTable) {

              if (theFinish == perfect) {
                  Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	          if (Rindex) {
                     Rindex2 = Rindex->Value(thePhotonMomentum);
                  }
                  else {
		     theStatus = NoRINDEX;
                     if ( verboseLevel > 0) BoundaryProcessVerbose();
                     aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                     aParticleChange.ProposeTrackStatus(fStopAndKill);
                     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
                  }
              }

              PropertyPointer =
                      aMaterialPropertiesTable->GetProperty("REFLECTIVITY");
              PropertyPointer1 =
                      aMaterialPropertiesTable->GetProperty("REALRINDEX");
              PropertyPointer2 =
                      aMaterialPropertiesTable->GetProperty("IMAGINARYRINDEX");

              iTE = 1;
              iTM = 1;

              if (PropertyPointer) {

                 theReflectivity =
                          PropertyPointer->Value(thePhotonMomentum);

              } else if (PropertyPointer1 && PropertyPointer2) {

                 CalculateReflectivity();

              }

              PropertyPointer =
              aMaterialPropertiesTable->GetProperty("EFFICIENCY");
              if (PropertyPointer) {
                      theEfficiency =
                      PropertyPointer->Value(thePhotonMomentum);
              }

              PropertyPointer =
              aMaterialPropertiesTable->GetProperty("TRANSMITTANCE");
              if (PropertyPointer) {
                      theTransmittance =
                      PropertyPointer->Value(thePhotonMomentum);
              }

              if (aMaterialPropertiesTable->
                                     ConstPropertyExists("SURFACEROUGHNESS"))
                 theSurfaceRoughness = aMaterialPropertiesTable->
                                         GetConstProperty("SURFACEROUGHNESS");

	      }
	   }
        else {
          // if not OpticalSurface}
          return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
        }

        if (type == dielectric_dielectric ) {
           if (theFinish == perfect ) {
              aMaterialPropertiesTable =
                     Material2->GetMaterialPropertiesTable();
              if (aMaterialPropertiesTable) {
                 Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
              }
              if (Rindex) {
                 Rindex2 = Rindex->Value(thePhotonMomentum);
              }
              else {
		 theStatus = NoRINDEX;
                 if ( verboseLevel > 0) BoundaryProcessVerbose();
                 aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                 aParticleChange.ProposeTrackStatus(fStopAndKill);
                 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	      }
           }
        }

	if (theModel == prism && theFinish == perfect) {
           G4double rand = G4UniformRand();
           if (rand > theReflectivity) {
              DoAbsorption();
           }
           else {
              DoRetroReflection();
           }
	}
	else {
	  G4cerr << " Error: WCSimOpRetroBoundaryProcess: illegal boundary type " << G4endl;
	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        NewMomentum = NewMomentum.unit();
        NewPolarization = NewPolarization.unit();

        if ( verboseLevel > 0) {
	   G4cout << " New Momentum Direction: " << NewMomentum     << G4endl;
	   G4cout << " New Polarization:       " << NewPolarization << G4endl;
           BoundaryProcessVerbose();
        }

	aParticleChange.ProposeMomentumDirection(NewMomentum);
	aParticleChange.ProposePolarization(NewPolarization);

        if ( theStatus == Detection ) InvokeSD(pStep);

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

void WCSimOpRetroBoundaryProcess::BoundaryProcessVerbose() const
{
        if ( theStatus == Undefined )
                G4cout << " *** Undefined *** " << G4endl;
        if ( theStatus == Transmission )
                G4cout << " *** Transmission *** " << G4endl;
        if ( theStatus == FresnelRefraction )
                G4cout << " *** FresnelRefraction *** " << G4endl;
        if ( theStatus == FresnelReflection )
                G4cout << " *** FresnelReflection *** " << G4endl;
        if ( theStatus == TotalInternalReflection )
                G4cout << " *** TotalInternalReflection *** " << G4endl;
        if ( theStatus == LambertianReflection )
                G4cout << " *** LambertianReflection *** " << G4endl;
        if ( theStatus == LobeReflection )
                G4cout << " *** LobeReflection *** " << G4endl;
        if ( theStatus == SpikeReflection )
                G4cout << " *** SpikeReflection *** " << G4endl;
        if ( theStatus == BackScattering )
                G4cout << " *** BackScattering *** " << G4endl;
        if ( theStatus == Absorption )
                G4cout << " *** Absorption *** " << G4endl;
        if ( theStatus == Detection )
                G4cout << " *** Detection *** " << G4endl;
        if ( theStatus == NotAtBoundary )
                G4cout << " *** NotAtBoundary *** " << G4endl;
        if ( theStatus == SameMaterial )
                G4cout << " *** SameMaterial *** " << G4endl;
        if ( theStatus == StepTooSmall )
                G4cout << " *** StepTooSmall *** " << G4endl;
        if ( theStatus == NoRINDEX )
                G4cout << " *** NoRINDEX *** " << G4endl;
}

G4ThreeVector
WCSimOpRetroBoundaryProcess::GetFacetNormal(const G4ThreeVector& Momentum,
			            const G4ThreeVector&  Normal ) const
{
        G4ThreeVector FacetNormal;

	   G4double  polish = 1.0;
	   if (OpticalSurface) polish = OpticalSurface->GetPolish();

           if (polish < 1.0) {
              do {
                 G4ThreeVector smear;
                 do {
                    smear.setX(2.*G4UniformRand()-1.0);
                    smear.setY(2.*G4UniformRand()-1.0);
                    smear.setZ(2.*G4UniformRand()-1.0);
                 } while (smear.mag()>1.0);
                 smear = (1.-polish) * smear;
                 FacetNormal = Normal + smear;
              } while (Momentum * FacetNormal >= 0.0);
              FacetNormal = FacetNormal.unit();
	   }
           else {
              FacetNormal = Normal;
           }
	return FacetNormal;
}


// GetMeanFreePath
// ---------------
//
G4double WCSimOpRetroBoundaryProcess::GetMeanFreePath(const G4Track& ,
                                              G4double ,
                                              G4ForceCondition* condition)
{
	*condition = Forced;

	return DBL_MAX;
}

G4double WCSimOpRetroBoundaryProcess::GetIncidentAngle() 
{
    G4double PdotN = OldMomentum * theFacetNormal;
    G4double magP= OldMomentum.mag();
    G4double magN= theFacetNormal.mag();
    G4double incidentangle = pi - std::acos(PdotN/(magP*magN));

    return incidentangle;
}

G4double WCSimOpRetroBoundaryProcess::GetReflectivity(G4double E1_perp,
                                              G4double E1_parl,
                                              G4double incidentangle,
                                              G4double RealRindex,
                                              G4double ImaginaryRindex)
{

  G4complex Reflectivity, Reflectivity_TE, Reflectivity_TM;
  G4complex N(RealRindex, ImaginaryRindex);
  G4complex CosPhi;

  G4complex u(1,0);           //unit number 1

  G4complex numeratorTE;      // E1_perp=1 E1_parl=0 -> TE polarization
  G4complex numeratorTM;      // E1_parl=1 E1_perp=0 -> TM polarization
  G4complex denominatorTE, denominatorTM;
  G4complex rTM, rTE;

  // Following two equations, rTM and rTE, are from: "Introduction To Modern
  // Optics" written by Fowles

  CosPhi=std::sqrt(u-((std::sin(incidentangle)*std::sin(incidentangle))/(N*N)));

  numeratorTE   = std::cos(incidentangle) - N*CosPhi;
  denominatorTE = std::cos(incidentangle) + N*CosPhi;
  rTE = numeratorTE/denominatorTE;

  numeratorTM   = N*std::cos(incidentangle) - CosPhi;
  denominatorTM = N*std::cos(incidentangle) + CosPhi;
  rTM = numeratorTM/denominatorTM;

  // This is my calculaton for reflectivity on a metalic surface
  // depending on the fraction of TE and TM polarization
  // when TE polarization, E1_parl=0 and E1_perp=1, R=abs(rTE)^2 and
  // when TM polarization, E1_parl=1 and E1_perp=0, R=abs(rTM)^2

  Reflectivity_TE =  (rTE*conj(rTE))*(E1_perp*E1_perp)
                    / (E1_perp*E1_perp + E1_parl*E1_parl);
  Reflectivity_TM =  (rTM*conj(rTM))*(E1_parl*E1_parl)
                    / (E1_perp*E1_perp + E1_parl*E1_parl);
  Reflectivity    = Reflectivity_TE + Reflectivity_TM;

  do {
     if(G4UniformRand()*real(Reflectivity) > real(Reflectivity_TE))
       {iTE = -1;}else{iTE = 1;}
     if(G4UniformRand()*real(Reflectivity) > real(Reflectivity_TM))
       {iTM = -1;}else{iTM = 1;}
  } while(iTE<0&&iTM<0);

  return real(Reflectivity);

}

void WCSimOpRetroBoundaryProcess::CalculateReflectivity()
{
  G4double RealRindex =
           PropertyPointer1->Value(thePhotonMomentum);
  G4double ImaginaryRindex =
           PropertyPointer2->Value(thePhotonMomentum);

  // calculate FacetNormal
  theFacetNormal = theGlobalNormal;

  G4double PdotN = OldMomentum * theFacetNormal;
  cost1 = -PdotN;

  if (std::abs(cost1) < 1.0 - kCarTolerance) {
     sint1 = std::sqrt(1. - cost1*cost1);
  } else {
     sint1 = 0.0;
  }

  G4ThreeVector A_trans, A_paral, E1pp, E1pl;
  G4double E1_perp, E1_parl;

  if (sint1 > 0.0 ) {
     A_trans = OldMomentum.cross(theFacetNormal);
     A_trans = A_trans.unit();
     E1_perp = OldPolarization * A_trans;
     E1pp    = E1_perp * A_trans;
     E1pl    = OldPolarization - E1pp;
     E1_parl = E1pl.mag();
  }
  else {
     A_trans  = OldPolarization;
     // Here we Follow Jackson's conventions and we set the
     // parallel component = 1 in case of a ray perpendicular
     // to the surface
     E1_perp  = 0.0;
     E1_parl  = 1.0;
  }

  //calculate incident angle
  G4double incidentangle = GetIncidentAngle();

  //calculate the reflectivity depending on incident angle,
  //polarization and complex refractive

  theReflectivity =
             GetReflectivity(E1_perp, E1_parl, incidentangle,
                                                 RealRindex, ImaginaryRindex);
}

G4bool WCSimOpRetroBoundaryProcess::InvokeSD(const G4Step* pStep)
{
  G4Step aStep = *pStep;

  aStep.AddTotalEnergyDeposit(thePhotonMomentum);

  G4VSensitiveDetector* sd = aStep.GetPostStepPoint()->GetSensitiveDetector();
  if (sd) return sd->Hit(&aStep);
  else return false;
}
