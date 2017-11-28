// 
// Modified from Geant4's source/processes/optical/include/WCSimOpRetroBoundaryProcess.hh
//
// Author:      Lukas Berns
// mail:        lukas.berns@hep.phys.titech.ac.jp
//
////////////////////////////////////////////////////////////////////////

#ifndef WCSimOpRetroBoundaryProcess_h
#define WCSimOpRetroBoundaryProcess_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"

#include "G4RandomTools.hh"
#include "G4RandomDirection.hh"

#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "WCSimOpticalSurface.hh"

#include "G4OpticalPhoton.hh"
#include "G4TransportationManager.hh"

// Class Description:
// Discrete Process -- retro-reflection at special optical interfaces.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

enum WCSimOpRetroBoundaryProcessStatus {
                                  Undefined,
                                  Transmission, FresnelRefraction,
                                  FresnelReflection, TotalInternalReflection,
                                  LambertianReflection, LobeReflection,
                                  SpikeReflection, BackScattering,
                                  RetroReflection, // new
                                  Absorption, Detection, NotAtBoundary,
                                  SameMaterial, StepTooSmall, NoRINDEX };

class WCSimOpRetroBoundaryProcess : public G4VDiscreteProcess
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        WCSimOpRetroBoundaryProcess(const G4String& processName = "OpRetroBoundary",
                                     G4ProcessType type = fOptical);
	    ~WCSimOpRetroBoundaryProcess();

private:

        WCSimOpRetroBoundaryProcess(const WCSimOpRetroBoundaryProcess &right);

        //////////////
        // Operators
        //////////////

        WCSimOpRetroBoundaryProcess& operator=(const WCSimOpRetroBoundaryProcess &right);

public:

	////////////
	// Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an optical photon.

	G4double GetMeanFreePath(const G4Track& ,
				 G4double ,
				 G4ForceCondition* condition);
        // Returns infinity; i. e. the process does not limit the step,
        // but sets the 'Forced' condition for the DoIt to be invoked at
        // every step. However, only at a boundary will any action be
        // taken.

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				       const G4Step&  aStep);
        // This is the method implementing boundary processes.

        WCSimOpRetroBoundaryProcessStatus GetStatus() const;

private:

	G4bool G4BooleanRand(const G4double prob) const;
    G4ThreeVector GetFacetNormal(const G4ThreeVector& Momentum,
                         const G4ThreeVector&  Normal) const;


        void ChooseReflection();
        void DoAbsorption();
        void DoReflection();
        void DoRetroReflection();

        G4double GetIncidentAngle();
        // Returns the incident angle of optical photon

        G4double GetReflectivity(G4double E1_perp,
                                 G4double E1_parl,
                                 G4double incidentangle,
                                 G4double RealRindex,
                                 G4double ImaginaryRindex);
        // Returns the Reflectivity on a metalic surface

        void CalculateReflectivity(void);

        void BoundaryProcessVerbose(void) const;

        // Invoke SD for post step point if the photon is 'detected'
        G4bool InvokeSD(const G4Step* step);

private:

	G4double thePhotonMomentum;

	G4ThreeVector OldMomentum;
	G4ThreeVector OldPolarization;

	G4ThreeVector NewMomentum;
	G4ThreeVector NewPolarization;

	G4ThreeVector theGlobalNormal;
	G4ThreeVector theFacetNormal;

	G4Material* Material1;
	G4Material* Material2;

	WCSimOpticalSurface* OpticalSurface;

        G4MaterialPropertyVector* PropertyPointer;
        G4MaterialPropertyVector* PropertyPointer1;
        G4MaterialPropertyVector* PropertyPointer2;

	G4double Rindex1;
	G4double Rindex2;

	G4double cost1, cost2, sint1, sint2;

	WCSimOpRetroBoundaryProcessStatus theStatus;

	WCSimOpticalSurfaceModel theModel;

	WCSimOpticalSurfaceFinish theFinish;

	G4double theReflectivity;
	G4double theEfficiency;
        G4double theTransmittance;

        G4double theSurfaceRoughness;

	G4double prob_sl, prob_ss, prob_bs;

        G4int iTE, iTM;

        G4double kCarTolerance;

        size_t idx, idy;
        G4Physics2DVector* DichroicVector;
};

////////////////////
// Inline methods
////////////////////

inline
G4bool WCSimOpRetroBoundaryProcess::G4BooleanRand(const G4double prob) const
{
  /* Returns a random boolean variable with the specified probability */

  return (G4UniformRand() < prob);
}

inline
G4bool WCSimOpRetroBoundaryProcess::IsApplicable(const G4ParticleDefinition& 
					               aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
void WCSimOpRetroBoundaryProcess::DoAbsorption()
{
              theStatus = Absorption;

              if ( G4BooleanRand(theEfficiency) ) {
		
                 // EnergyDeposited =/= 0 means: photon has been detected
                 theStatus = Detection;
                 aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
              }
              else {
                 aParticleChange.ProposeLocalEnergyDeposit(0.0);
              }

              NewMomentum = OldMomentum;
              NewPolarization = OldPolarization;

//              aParticleChange.ProposeEnergy(0.0);
              aParticleChange.ProposeTrackStatus(fStopAndKill);
}

inline
void WCSimOpRetroBoundaryProcess::DoReflection()
{

        theStatus = SpikeReflection;
          theFacetNormal = theGlobalNormal;
          G4double PdotN = OldMomentum * theFacetNormal;
          NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

        G4double EdotN = OldPolarization * theFacetNormal;
        NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
}

inline
void WCSimOpRetroBoundaryProcess::DoRetroReflection()
{
    theStatus = RetroReflection;
    NewMomentum = -OldMomentum;
    NewPolarization = OldPolarization;
}

#endif /* WCSimOpRetroBoundaryProcess_h */
