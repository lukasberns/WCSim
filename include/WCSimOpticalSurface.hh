// Modified from Geant4's WCSimOpticalSurface
// Author:      Lukas Berns
// mail:        lukas.berns@hep.phys.titech.ac.jp
//
////////////////////////////////////////////////////////////////////////

#ifndef WCSimOpticalSurface_h
#define WCSimOpticalSurface_h 1

/////////////
// Includes
/////////////

#include "G4Types.hh"
#include "G4Physics2DVector.hh"
#include "G4SurfaceProperty.hh"

// Class Description:
// A optical surface class for use in the WCSimOpRetroBoundaryProcess class.
// Contains the enumerations: WCSimOpticalSurfaceFinish, WCSimOpticalSurfaceType,
// and WCSimOpticalSurfaceModel.
// Class Description - End:

enum WCSimOpticalSurfaceFinish
{
   perfect
};

enum WCSimOpticalSurfaceModel
{
   prism
};

class G4MaterialPropertiesTable;

/////////////////////
// Class Definition
/////////////////////

class WCSimOpticalSurface : public G4SurfaceProperty
{

public: // Without description
  
        //////////////
        // Operators
        //////////////
  
	WCSimOpticalSurface(const WCSimOpticalSurface &right);
	WCSimOpticalSurface & operator=(const WCSimOpticalSurface &right);
  
	G4int operator==(const WCSimOpticalSurface &right) const;
	G4int operator!=(const WCSimOpticalSurface &right) const;

public: // With description

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

	WCSimOpticalSurface(const G4String& name,
			 WCSimOpticalSurfaceModel model = prism,
			 WCSimOpticalSurfaceFinish finish = perfect,
			 G4SurfaceType type = dielectric_dielectric,
			 G4double value = 1.0);
        // Constructor of an optical surface object.

public: // Without description

	virtual ~WCSimOpticalSurface();

	////////////
	// Methods
        ////////////

	// public methods

public: // With description

        void SetType(const G4SurfaceType& type);

        inline WCSimOpticalSurfaceFinish GetFinish() const { return theFinish; }
        // Returns the optical surface finish.
        void SetFinish(const WCSimOpticalSurfaceFinish );
        // Sets the optical surface finish.

        inline WCSimOpticalSurfaceModel GetModel() const { return theModel; }
        // Returns the optical surface model used.
        inline void SetModel(const WCSimOpticalSurfaceModel model)
                                                      { theModel = model; }
        // Sets the optical surface model to be followed.

	inline G4double GetSigmaAlpha() const { return sigma_alpha; }
        // Returns an unified model surface parameter.
	inline void SetSigmaAlpha(const G4double s_a) { sigma_alpha = s_a; }
        // Sets an unified model surface parameter.

	G4double GetPolish() const { return polish; }
        // Returns the optical surface polish type.
	inline void SetPolish(const G4double plsh) { polish=plsh; }
        // Sets the optical surface polish type.

	inline G4MaterialPropertiesTable* GetMaterialPropertiesTable() const
				       { return theMaterialPropertiesTable; }
        // Retrieves the pointer of the G4MaterialPropertiesTable 
        // attached to optical surface.

	inline void SetMaterialPropertiesTable(G4MaterialPropertiesTable *anMPT)
				       { theMaterialPropertiesTable = anMPT; }
        // Attaches a G4MaterialPropertiesTable to the optical surface.

	void DumpInfo() const;
        // Prints information about the optical surface.

        void ReadLUTFile(void);
        // Method to read the Look-Up-Table into array AngularDistribution

        inline G4double GetAngularDistributionValue(G4int, G4int, G4int);

        inline G4int GetThetaIndexMax(void) const { return thetaIndexMax; }
        inline G4int GetPhiIndexMax(void) const { return phiIndexMax; }

        void ReadDichroicFile(void);
        // Method to read the dichroic surface data file into Dichroic

        inline G4Physics2DVector* GetDichroicVector();

private:

// ------------------
// Basic data members ( To define an optical surface)
// ------------------

        WCSimOpticalSurfaceModel theModel;		// Surface model
        WCSimOpticalSurfaceFinish theFinish;	// Surface finish

	G4double sigma_alpha;		// The sigma of micro-facet polar angle
	G4double polish;		// Polish parameter in glisur model

	G4MaterialPropertiesTable* theMaterialPropertiesTable;

        static const G4int incidentIndexMax = 91;
        static const G4int thetaIndexMax = 45;
        static const G4int phiIndexMax = 37;

        G4float* AngularDistribution;

        G4Physics2DVector* DichroicVector;

};

////////////////////
// Inline methods
////////////////////

inline
 G4double WCSimOpticalSurface::GetAngularDistributionValue(G4int angleIncident,
                                                        G4int thetaIndex,
                                                        G4int phiIndex)
{
  return AngularDistribution[angleIncident+
                             thetaIndex*incidentIndexMax+
                             phiIndex*thetaIndexMax*incidentIndexMax];
}

inline
 G4Physics2DVector* WCSimOpticalSurface::GetDichroicVector()
{
  return DichroicVector;
}

#endif /* WCSimOpticalSurface_h */
