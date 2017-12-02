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

#include "WCSimOpticalRetroSurface.hh"

#include "G4OpticalPhoton.hh"
#include "G4TransportationManager.hh"
#include "G4OpBoundaryProcess.hh"

// Class Description:
// Discrete Process -- retro-reflection at special optical interfaces.
// Class inherits publicly from G4OpBoundaryProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

enum WCSimOpRetroBoundaryProcessStatus {
                                  RetroUndefined,
                                  RetroReflection
                                  };

class WCSimOpRetroBoundaryProcess : public G4OpBoundaryProcess
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
        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                        const G4Step&  aStep);


        WCSimOpRetroBoundaryProcessStatus GetRetroStatus() const;

protected:

        void DoRetroReflection();

        // overridden method
        virtual void CustomBoundary();
        void BoundaryProcessVerbose(void) const;

protected:

	WCSimOpticalRetroSurface* OpticalRetroSurface;
	WCSimOpRetroBoundaryProcessStatus theRetroStatus;

        // access to superclass variables
        
};

////////////////////
// Inline methods
////////////////////

inline
void WCSimOpRetroBoundaryProcess::DoRetroReflection()
{
    theRetroStatus = RetroReflection;
    NewMomentum = -OldMomentum;
    NewPolarization = OldPolarization;
}

#endif /* WCSimOpRetroBoundaryProcess_h */
