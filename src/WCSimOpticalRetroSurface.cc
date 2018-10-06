// Modified from Geant4's G4OpticalSurface
// Author:      Lukas Berns
// mail:        lukas.berns@hep.phys.titech.ac.jp
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//#include "G4ios.hh"
#include "globals.hh"
#include "WCSimOpticalRetroSurface.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

WCSimOpticalRetroSurface& WCSimOpticalRetroSurface::operator=(const WCSimOpticalRetroSurface& right)
{
  G4OpticalSurface::operator=(right);
  if (this != &right)
    {
      theRetroModel                   = right.theRetroModel;
     } 
  return *this;
}

        /////////////////
        // Constructors
        /////////////////

WCSimOpticalRetroSurface::WCSimOpticalRetroSurface(const G4String& name,
				   WCSimOpticalRetroSurfaceModel model,
				   G4double value)
                                   : G4OpticalSurface(name,glisur,polished,custom_boundary),
				     theRetroModel(model)
{
	if (model == prism ){
	}
	else {
                G4Exception("WCSimOpticalSurface::WCSimOpticalSurface()", "mat309",
                            FatalException,
			    "Constructor called with INVALID model.");
	}
}

WCSimOpticalRetroSurface::~WCSimOpticalRetroSurface()
{
}

WCSimOpticalRetroSurface::WCSimOpticalRetroSurface(const WCSimOpticalRetroSurface &right)
  : G4OpticalSurface(right)
{
       *this = right;
       this->theRetroModel = right.theRetroModel;
}

G4int WCSimOpticalRetroSurface::operator==(const WCSimOpticalRetroSurface &right) const
{
	return (this == (WCSimOpticalRetroSurface *) &right);
}

G4int WCSimOpticalRetroSurface::operator!=(const WCSimOpticalRetroSurface &right) const
{
	return (this != (WCSimOpticalRetroSurface *) &right);
}
        ////////////
        // Methods
        ////////////

