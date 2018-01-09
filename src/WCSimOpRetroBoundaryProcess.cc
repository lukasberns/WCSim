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
        return G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep);
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
                if (rand > theReflectivity*theReflectivity*theReflectivity) { // three reflections, so ^3
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

