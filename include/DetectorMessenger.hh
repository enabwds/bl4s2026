#pragma once

// ============================================================
//  DetectorMessenger.hh
//
//  UI commands for configuring the absorber.
//
//  Available commands:
//
//    /det/setAbsorberMaterial <material>
//    /det/setAbsorberThickness <value> <unit>
//    /det/update
//
//  The messenger only passes configuration changes to
//  DetectorConstruction. It does not contain physics or
//  analysis logic.
// ============================================================

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class DetectorMessenger : public G4UImessenger
{
public:

    explicit DetectorMessenger(DetectorConstruction* detector);

    ~DetectorMessenger() override;

    void SetNewValue(G4UIcommand* command,
                     G4String newValue) override;

private:

    DetectorConstruction* fDetector = nullptr;

    G4UIcmdWithAString*        fMaterialCommand   = nullptr;
    G4UIcmdWithADoubleAndUnit* fThicknessCommand  = nullptr;
    G4UIcmdWithoutParameter*   fUpdateCommand     = nullptr;
};

