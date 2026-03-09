// ============================================================
//  DetectorMessenger.cc
// ============================================================

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
: fDetector(det)
{
    fDetDir = new G4UIdirectory("/det/");
    fDetDir->SetGuidance("Detector control commands.");

    fMatCmd = new G4UIcmdWithAString("/det/setAbsorberMaterial", this);
    fMatCmd->SetGuidance("Set the absorber material (NIST name, e.g. G4_Fe).");
    fMatCmd->SetParameterName("Material", false);
    fMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fThickCmd = new G4UIcmdWithADoubleAndUnit("/det/setAbsorberThickness", this);
    fThickCmd->SetGuidance("Set absorber slab thickness.");
    fThickCmd->SetParameterName("Thickness", false);
    fThickCmd->SetUnitCategory("Length");
    fThickCmd->SetDefaultUnit("mm");
    fThickCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fUpdateCmd = new G4UIcmdWithoutParameter("/det/update", this);
    fUpdateCmd->SetGuidance("Update geometry after parameter changes.");
    fUpdateCmd->AvailableForStates(G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
    delete fMatCmd;
    delete fThickCmd;
    delete fUpdateCmd;
    delete fDetDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
    if (command == fMatCmd)
        fDetector->SetAbsorberMaterial(value);

    else if (command == fThickCmd)
        fDetector->SetAbsorberThickness(
            fThickCmd->GetNewDoubleValue(value));

    else if (command == fUpdateCmd)
        G4RunManager::GetRunManager()->ReinitializeGeometry();
}
