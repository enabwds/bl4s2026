#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // --------------------------------------------------------
    // Detector configuration
    // --------------------------------------------------------

    void SetAbsorberMaterial(const G4String& materialName);
    void SetAbsorberThickness(G4double thickness);

    void UpdateGeometry();

    // --------------------------------------------------------
    // Detector information
    // --------------------------------------------------------

    G4Material* GetAbsorberMaterial() const;
    G4double GetAbsorberThickness() const;

private:
    void DefineMaterials();
    void PrintParameters() const;

    // --------------------------------------------------------
    // Geometry constants
    // --------------------------------------------------------

    static constexpr G4int kNCols = 4;
    static constexpr G4int kNRows = 4;

    static constexpr G4double kBlockSizeXY = 10.0 * cm;
    static constexpr G4double kBlockSizeZ  = 37.0 * cm;

    static constexpr G4double kAbsorberSizeXY = 20.0 * cm;

    // Downstream face of absorber
    static constexpr G4double kAbsorberDownstreamZ = -6.0 * cm;

    // --------------------------------------------------------
    // Mutable detector configuration
    // --------------------------------------------------------

    G4Material* fAbsorberMaterial = nullptr;

    // Store with GEANT4 units, not as a bare number.
    G4double fAbsorberThickness = 17.57 * mm;

    // --------------------------------------------------------
    // Geometry
    // --------------------------------------------------------

    G4LogicalVolume* fAbsorberLogical = nullptr;
    G4LogicalVolume* fScoringVolume = nullptr;

    DetectorMessenger* fMessenger = nullptr;
};

#endif