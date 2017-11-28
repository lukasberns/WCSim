// Modified from Geant4's G4OpticalSurface
// Author:      Lukas Berns
// mail:        lukas.berns@hep.phys.titech.ac.jp
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//#include "G4ios.hh"
#include "globals.hh"
#include "WCSimOpticalSurface.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

WCSimOpticalSurface& WCSimOpticalSurface::operator=(const WCSimOpticalSurface& right)
{
  if (this != &right)
    {
      theName                    = right.theName;
      theType                    = right.theType;
      theModel                   = right.theModel;
      theFinish                  = right.theFinish;
      sigma_alpha                = right.sigma_alpha;
      polish                     = right.polish;
      theMaterialPropertiesTable = right.theMaterialPropertiesTable;
      if (AngularDistribution) delete [] AngularDistribution;
      AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
      *(AngularDistribution)     = *(right.AngularDistribution);
      if (DichroicVector) delete DichroicVector;
      DichroicVector = new G4Physics2DVector();
      *DichroicVector            = *(right.DichroicVector);
     } 
  return *this;
}

        /////////////////
        // Constructors
        /////////////////

WCSimOpticalSurface::WCSimOpticalSurface(const G4String& name,
				   WCSimOpticalSurfaceModel model,
				   WCSimOpticalSurfaceFinish finish,
				   G4SurfaceType type,
				   G4double value)
                                   : G4SurfaceProperty(name,type),
				     theModel(model),
				     theFinish(finish),
				     theMaterialPropertiesTable(0)
{
	if (model == prism ){
		polish = value;
		sigma_alpha = 0.0;
	}
	else {
                G4Exception("WCSimOpticalSurface::WCSimOpticalSurface()", "mat309",
                            FatalException,
			    "Constructor called with INVALID model.");
	}

        AngularDistribution = NULL;
        DichroicVector      = NULL;

}

WCSimOpticalSurface::~WCSimOpticalSurface()
{
        if (AngularDistribution) delete [] AngularDistribution;
        if (DichroicVector) delete DichroicVector;
}

WCSimOpticalSurface::WCSimOpticalSurface(const WCSimOpticalSurface &right)
  : G4SurfaceProperty(right.theName,right.theType)
{
       *this = right;
       this->theName = right.theName;
       this->theType = right.theType;
       this->theModel = right.theModel;
       this->theFinish = right.theFinish;
       this->sigma_alpha = right.sigma_alpha;
       this->polish = right.polish;
       this->theMaterialPropertiesTable = right.theMaterialPropertiesTable;
       if (this->AngularDistribution) delete [] AngularDistribution;
       this->AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
       *(this->AngularDistribution) = *(right.AngularDistribution);
       if (this->DichroicVector) delete DichroicVector;
       this->DichroicVector = new G4Physics2DVector();
       *(this->DichroicVector) = *(right.DichroicVector);
}

G4int WCSimOpticalSurface::operator==(const WCSimOpticalSurface &right) const
{
	return (this == (WCSimOpticalSurface *) &right);
}

G4int WCSimOpticalSurface::operator!=(const WCSimOpticalSurface &right) const
{
	return (this != (WCSimOpticalSurface *) &right);
}
        ////////////
        // Methods
        ////////////

void WCSimOpticalSurface::DumpInfo() const 
{

	// Dump info for surface

	G4cout << 
        "  Surface type   = " << G4int(theType)   << G4endl <<
	"  Surface finish = " << G4int(theFinish) << G4endl <<
	"  Surface model  = " << G4int(theModel)  << G4endl;

	G4cout << G4endl;
}

void WCSimOpticalSurface::SetType(const G4SurfaceType& type)
{
  theType = type;
}

void WCSimOpticalSurface::SetFinish(const WCSimOpticalSurfaceFinish finish)
{
  theFinish = finish;
}

void WCSimOpticalSurface::ReadLUTFile()
{
}

void WCSimOpticalSurface::ReadDichroicFile()
{
  const char* datadir = getenv("G4DICHROICDATA");

  if(!datadir) {
    G4Exception("WCSimOpticalSurface::ReadDichroicFile()","mat313",
    FatalException,"Environment variable G4DICHROICDATA not defined");
    return;
  }

  std::ostringstream ost;
  ost << datadir;
  std::ifstream fin(ost.str().c_str());
  if( !fin.is_open()) {
    G4ExceptionDescription ed;
    ed << "Dichroic surface data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("WCSimOpticalSurface::ReadDichroicFile()","mat314",
    FatalException,ed," ");
    return;
  }

  if( !(DichroicVector->Retrieve(fin)) ) {
    G4ExceptionDescription ed;
    ed << "Dichroic surface data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("WCSimOpticalSurface::ReadDichroicFile()","mat315",
    FatalException,ed," ");
    return;
  }

//  DichroicVector->SetBicubicInterpolation(true);

  G4cout << " *** Dichroic surface data file *** " << G4endl;

  G4int numberOfXNodes = DichroicVector->GetLengthX();
  G4int numberOfYNodes = DichroicVector->GetLengthY();

  G4cout << "numberOfXNodes: " << numberOfXNodes << G4endl;
  G4cout << "numberOfYNodes: " << numberOfYNodes << G4endl;

  if (0 > numberOfXNodes || numberOfXNodes >= INT_MAX) numberOfXNodes = 0;
  if (0 > numberOfYNodes || numberOfYNodes >= INT_MAX) numberOfYNodes = 0;

  G4PV2DDataVector  xVector;
  G4PV2DDataVector  yVector;

  xVector.resize(numberOfXNodes,0.);
  yVector.resize(numberOfYNodes,0.);

  for(G4int i = 0; i<numberOfXNodes; ++i) {
     G4cout << "i: " << DichroicVector->GetX(i) << G4endl;
     xVector[i] = DichroicVector->GetX(i);
  }
  for(G4int j = 0; j<numberOfYNodes; ++j) {
     G4cout << "j: " << DichroicVector->GetY(j) << G4endl;
     yVector[j] = DichroicVector->GetY(j);
  }

  for(G4int j = 0; j<numberOfYNodes; ++j) {
     for(G4int i = 0; i<numberOfXNodes; ++i) {
        G4cout << " i: " << i << " j: " << j << " "
               << DichroicVector->GetValue(i,j) << G4endl;
     }
  }

//  G4int idx, idy;

//  for(G4int j = 0; j<numberOfYNodes-1; ++j) {
//     G4double y = (yVector[j] + yVector[j+1])/2.;
//     for(G4int i = 0; i<numberOfXNodes-1; ++i) {
//        G4double x = (xVector[i] + xVector[i+1])/2.;
//        G4cout << " x: " << x << " y: " << y << " "
//               << DichroicVector->Value(x,y,idx,idy) << G4endl;
//     }
//  }

}
