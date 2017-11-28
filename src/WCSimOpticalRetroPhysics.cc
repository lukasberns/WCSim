// Modified from Geant4's G4OpticalPhysics
//
// ClassName:   WCSimOpticalRetroPhysics
//
// Author:      Lukas Berns 27.11.2017
//
//----------------------------------------------------------------------------
//

#include "WCSimOpticalRetroPhysics.hh"

#include "WCSimOpRetroBoundaryProcess.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(WCSimOpticalRetroPhysics);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCSimOpticalRetroPhysics::WCSimOpticalRetroPhysics(G4int verbose, const G4String& name)
  : G4VPhysicsConstructor(name)
{
  verboseLevel = verbose;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCSimOpticalRetroPhysics::~WCSimOpticalRetroPhysics()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WCSimOpticalRetroPhysics::PrintStatistics() const
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WCSimOpticalRetroPhysics::ConstructParticle()
{
/// Instantiate particles.

  // the optical photon should be defined in G4OpticalPhysics
}

#include "G4Threading.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4AutoDelete.hh"

void WCSimOpticalRetroPhysics::ConstructProcess()
{
// Construct optical processes.

  if(verboseLevel>0)
         G4cout <<"WCSimOpticalRetroPhysics:: Add Optical Physics Processes"<< G4endl;

  // Add Optical Processes

  WCSimOpRetroBoundaryProcess* OpBoundaryProcess = new WCSimOpRetroBoundaryProcess();

  G4ProcessManager * pManager = 0;
  pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (!pManager) {
     std::ostringstream o;
     o << "Optical Photon without a Process Manager";
     G4Exception("WCSimOpticalRetroPhysics::ConstructProcess()","",
                  FatalException,o.str().c_str());
     return;
  }

  pManager->AddDiscreteProcess(OpBoundaryProcess);
  OpBoundaryProcess->SetVerboseLevel(verboseLevel);

  if (verboseLevel > 1) PrintStatistics();
  if (verboseLevel > 0)
    G4cout << "### " << namePhysics << " physics constructed." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
