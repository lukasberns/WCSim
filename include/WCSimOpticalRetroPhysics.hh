// Modified from Geant4's G4OpticalPhysics
// ClassName:   WCSimOpticalRetroPhysics
//
// Author:      Lukas Berns 27.11.2017
//
//---------------------------------------------------------------------------
//
// This class provides construction of retro-reflector optical physics
//

#ifndef WCSimOpticalRetroPhysics_h
#define WCSimOpticalRetroPhysics_h 1

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4OpticalProcessIndex.hh"

#include "WCSimOpticalSurface.hh"

#include <vector>

class G4VProcess;
class G4EmSaturation;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class WCSimOpticalRetroPhysics : public G4VPhysicsConstructor
{
  public:

    WCSimOpticalRetroPhysics(G4int verbose = 0, const G4String& name = "Optical");
    virtual ~WCSimOpticalRetroPhysics();

  protected:

    // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:

    /// Not implemented
    WCSimOpticalRetroPhysics(const WCSimOpticalRetroPhysics& right);
    /// Not implemented
    WCSimOpticalRetroPhysics& operator=(const WCSimOpticalRetroPhysics& right);

  public:

    // configure WCSimOpticalRetroPhysics builder
    void Configure(G4OpticalProcessIndex, G4bool );

  private:

    void PrintStatistics() const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // WCSimOpticalRetroPhysics_h
