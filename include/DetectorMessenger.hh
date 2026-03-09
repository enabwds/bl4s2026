#pragma once
// ============================================================
//  DetectorMessenger.hh
//  Exposes detector geometry parameters to GEANT4 macro commands.
//
//  Usage in .mac files:
//    /det/setAbsorberMaterial G4_Fe
//    /det/setAbsorberThickness 20 mm
//    /det/update
// ============================================================

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

class DetectorConstruction;

class DetectorMessenger : public G4UImessenger
{
public:
    explicit DetectorMessenger(DetectorConstruction* det);
    ~DetectorMessenger() override;

    void SetNewValue(G4UIcommand* command, G4String value) override;

private:
    DetectorConstruction*          fDetector;
    G4UIdirectory*                 fDetDir;
    G4UIcmdWithAString*            fMatCmd;
    G4UIcmdWithADoubleAndUnit*     fThickCmd;
    G4UIcmdWithoutParameter*       fUpdateCmd;
};
