// ============================================================
//  DetectorMessenger.cc
// ============================================================

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"


// ============================================================
//  Constructor
// ============================================================

DetectorMessenger::DetectorMessenger(
    DetectorConstruction* detector)
    : fDetector(detector)
{
    // --------------------------------------------------------
    //  /det/
    // --------------------------------------------------------

    auto* detectorDirectory = new G4UIdirectory("/det/");

    detectorDirectory->SetGuidance(
        "Detector configuration commands."
    );


    // --------------------------------------------------------
    //  /det/setAbsorberMaterial
    // --------------------------------------------------------

    fMaterialCommand =
        new G4UIcmdWithAString(
            "/det/setAbsorberMaterial",
            this
        );

    fMaterialCommand->SetGuidance(
        "Set the absorber material."
    );

    fMaterialCommand->SetParameterName(
        "material",
        false
    );

    fMaterialCommand->AvailableForStates(
        G4State_PreInit,
        G4State_Idle
    );


    // --------------------------------------------------------
    //  /det/setAbsorberThickness
    // --------------------------------------------------------

    fThicknessCommand =
        new G4UIcmdWithADoubleAndUnit(
            "/det/setAbsorberThickness",
            this
        );

    fThicknessCommand->SetGuidance(
        "Set the physical absorber thickness."
    );

    fThicknessCommand->SetParameterName(
        "thickness",
        false
    );

    fThicknessCommand->SetDefaultUnit("mm");

    fThicknessCommand->SetRange(
        "thickness > 0."
    );

    fThicknessCommand->AvailableForStates(
        G4State_PreInit,
        G4State_Idle
    );


    // --------------------------------------------------------
    //  /det/update
    // --------------------------------------------------------

    fUpdateCommand =
        new G4UIcmdWithoutParameter(
            "/det/update",
            this
        );

    fUpdateCommand->SetGuidance(
        "Apply the current detector configuration."
    );

    fUpdateCommand->AvailableForStates(
        G4State_Idle
    );
}


// ============================================================
//  Destructor
// ============================================================

DetectorMessenger::~DetectorMessenger()
{
    delete fMaterialCommand;
    delete fThicknessCommand;
    delete fUpdateCommand;
}


// ============================================================
//  Command handling
// ============================================================

void DetectorMessenger::SetNewValue(
    G4UIcommand* command,
    G4String newValue)
{
    if (command == fMaterialCommand)
    {
        fDetector->SetAbsorberMaterial(newValue);
    }
    else if (command == fThicknessCommand)
    {
        fDetector->SetAbsorberThickness(
            fThicknessCommand->GetNewDoubleValue(newValue)
        );
    }
    else if (command == fUpdateCommand)
    {
        fDetector->UpdateGeometry();
    }
}
