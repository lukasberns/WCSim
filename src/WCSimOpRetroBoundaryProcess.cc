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
             : G4OpBoundaryProcess(processName, type)
{
	theRetroStatus = RetroUndefined;
        OpticalRetroSurface = NULL;
        fTrack = NULL;
        fStep = NULL;
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
        theRetroStatus = RetroUndefined;
        fTrack = &aTrack;
        fStep = &aStep;
        G4VParticleChange *aChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep);
        fTrack = NULL;
        fStep = NULL;
        return aChange;
}

void WCSimOpRetroBoundaryProcess::BoundaryProcessVerbose() const
{
        if ( theRetroStatus == RetroReflection )
                G4cout << " *** RetroReflection *** " << G4endl;
        else
                G4OpBoundaryProcess::BoundaryProcessVerbose();
}

void WCSimOpRetroBoundaryProcess::CustomBoundary()
{
        // override
	OpticalRetroSurface = 
           dynamic_cast <WCSimOpticalRetroSurface*> (OpticalSurface);
        if (OpticalRetroSurface) {
                theRetroStatus = RetroReflection;
                G4double rand = G4UniformRand();
                double incidence = OldMomentum.angle(-theGlobalNormal);
                double acceptance = (60.*deg - incidence) / (60.*deg - 20.*deg);
                if (acceptance > 1.) { acceptance = 1.; }
                if (acceptance < 0.) { acceptance = 0.; }
                if (rand > acceptance*theReflectivity*theReflectivity*theReflectivity) { // three reflections, so ^3
                   DoAbsorption();
                }
                else {
                   DoRetroReflection();
                }
        }
        else {
            G4cout << "Error: was expecting WCSimOpticalRetroSurface for custom_boundary" << G4endl;
        }
}

void WCSimOpRetroBoundaryProcess::DoRetroReflection() {
    // we kill the old particle and create a new one, so this does not get identified as a "scattering" process
    
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    G4StepPoint* pPostStepPoint = fStep->GetPostStepPoint();

    NewMomentum = -OldMomentum; // only the direction
    NewPolarization = OldPolarization;

    // Generate a new photon:
    G4DynamicParticle* aRetroPhoton =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), NewMomentum);
    aRetroPhoton->SetKineticEnergy(thePhotonMomentum);
    aRetroPhoton->SetPolarization
           (NewPolarization.x(),
            NewPolarization.y(),
            NewPolarization.z());

    G4double aSecondaryTime = pPostStepPoint->GetGlobalTime();
    G4ThreeVector aSecondaryPosition = pPostStepPoint->GetPosition();

    G4Track* aSecondaryTrack =
      new G4Track(aRetroPhoton,aSecondaryTime,aSecondaryPosition);
    aSecondaryTrack->SetTouchableHandle(fTrack->GetTouchableHandle());
    // aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);

    aSecondaryTrack->SetParentID(fTrack->GetTrackID());
    aParticleChange.AddSecondary(aSecondaryTrack);
}

