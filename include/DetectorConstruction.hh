#pragma once
// ============================================================
//  DetectorConstruction.hh
//  Full BL4S beamline geometry:
//    • 2× Threshold Cherenkov detectors (electron ID)
//    • 2× Delay Wire Chambers (beam tracking)
//    • 2× Scintillator triggers (coincidence)
//    • Absorber slab (variable material + thickness)
//    • 4×4 lead-glass calorimeter (10×10×37 cm blocks)
// ============================================================

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    void SetAbsorberMaterial(const G4String& name);
    void SetAbsorberThickness(G4double thickness);  // must be in G4 internal units (pass value * mm)

    const G4Material* GetAbsorberMaterial() const { return fAbsorberMaterial; }
    // Returns thickness in G4 internal units. Divide by mm to get mm value.
    G4double          GetAbsorberThickness() const { return fAbsorberThickness; }

private:
    void DefineMaterials();
    void PrintParameters() const;

    // Absorber slab — fAbsorberThickness stored in G4 internal units (mm=1).
    // Default: 17.57 mm = 1 X0 of iron (X0_Fe = 1.757 cm).
    // The literal 17.57 is already in internal units because mm=1 in GEANT4,
    // so no unit header is needed here.
    G4Material* fAbsorberMaterial  = nullptr;
    G4double    fAbsorberThickness = 17.57;  // G4 internal units (= 17.57 mm)

    // Real BL4S calorimeter block dimensions
    static constexpr int    kNcols       = 4;
    static constexpr int    kNrows       = 4;
    static constexpr double kBlockSizeXY = 10.0;  // cm — real BL4S spec
    static constexpr double kBlockSizeZ  = 37.0;  // cm — ~18 X0 in lead glass
    // Downstream face of absorber anchored here — keeps air gap to calorimeter
    // constant regardless of thickness (without this it drifts ~4cm across the
    // thickness scan, introducing a systematic unrelated to absorber material).
    static constexpr double kAbsorberDownstreamZ = -6.0;  // cm

    G4LogicalVolume* fScoringVolume = nullptr;
    DetectorMessenger* fMessenger   = nullptr;
};
