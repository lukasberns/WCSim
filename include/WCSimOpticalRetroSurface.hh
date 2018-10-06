// Modified from Geant4's WCSimOpticalSurface
// Author:      Lukas Berns
// mail:        lukas.berns@hep.phys.titech.ac.jp
//
////////////////////////////////////////////////////////////////////////

#ifndef WCSimOpticalRetroSurface_h
#define WCSimOpticalRetroSurface_h 1

/////////////
// Includes
/////////////

#include "G4Types.hh"
#include "G4Physics2DVector.hh"
#include "G4SurfaceProperty.hh"
#include "G4OpticalSurface.hh"

// Class Description:
// A optical surface class for use in the WCSimOpRetroBoundaryProcess class.
// Contains the enumerations: WCSimOpticalRetroSurfaceModel.
// Class Description - End:

enum WCSimOpticalRetroSurfaceModel
{
   prism
};

/////////////////////
// Class Definition
/////////////////////

class WCSimOpticalRetroSurface : public G4OpticalSurface
{

public: // Without description
  
        //////////////
        // Operators
        //////////////
  
	WCSimOpticalRetroSurface(const WCSimOpticalRetroSurface &right);
	WCSimOpticalRetroSurface & operator=(const WCSimOpticalRetroSurface &right);
  
	G4int operator==(const WCSimOpticalRetroSurface &right) const;
	G4int operator!=(const WCSimOpticalRetroSurface &right) const;

public: // With description

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

	WCSimOpticalRetroSurface(const G4String& name,
			 WCSimOpticalRetroSurfaceModel model = prism,
			 G4double value = 1.0);
        // Constructor of an optical surface object.

public: // Without description

	virtual ~WCSimOpticalRetroSurface();

	////////////
	// Methods
        ////////////

	// public methods

public: // With description

        inline WCSimOpticalRetroSurfaceModel GetRetroModel() const { return theRetroModel; }
        // Returns the optical surface model used.
        inline void SetRetroModel(const WCSimOpticalRetroSurfaceModel model)
                                                      { theRetroModel = model; }
        // Sets the optical surface model to be followed.

private:

// ------------------
// Basic data members ( To define an optical surface)
// ------------------

        WCSimOpticalRetroSurfaceModel theRetroModel;		// Surface model
};

////////////////////
// Inline methods
////////////////////

#endif /* WCSimOpticalRetroSurface_h */
