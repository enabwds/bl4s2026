#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "CalorimeterSD.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4StateManager.hh"
#include "G4RunManager.hh"
#include "G4Exception.hh"
#include "G4LogicalVolumeStore.hh"

#include <iomanip>
#include <iostream>

// ============================================================
// Constructor / Destructor
// ============================================================

DetectorConstruction::DetectorConstruction()
{
    DefineMaterials();

    fAbsorberMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");

    fMessenger = new DetectorMessenger(this);
}


DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}


// ============================================================
// Material definitions
// ============================================================

void DetectorConstruction::DefineMaterials()
{
    auto* nist = G4NistManager::Instance();

    // --------------------------------------------------------
    // Standard materials
    // --------------------------------------------------------

    nist->FindOrBuildMaterial("G4_AIR");
    nist->FindOrBuildMaterial("G4_Al");
    nist->FindOrBuildMaterial("G4_Fe");
    nist->FindOrBuildMaterial("G4_Cu");
    nist->FindOrBuildMaterial("G4_Pb");
    nist->FindOrBuildMaterial("G4_W");
    nist->FindOrBuildMaterial("G4_Sn");
    nist->FindOrBuildMaterial("G4_Ti");

    // --------------------------------------------------------
    // Lead glass
    //
    // Composition used by the existing simulation.
    // --------------------------------------------------------

    if (!G4Material::GetMaterial("LeadGlass", false))
    {
        auto* Pb = nist->FindOrBuildMaterial("G4_Pb");
        auto* Si = nist->FindOrBuildMaterial("G4_Si");
        auto* O  = nist->FindOrBuildMaterial("G4_O");

        auto* leadGlass =
            new G4Material(
                "LeadGlass",
                3.86 * g / cm3,
                3
            );

        leadGlass->AddMaterial(Pb, 0.474);
        leadGlass->AddMaterial(Si, 0.229);
        leadGlass->AddMaterial(O,  0.297);
    }

    // --------------------------------------------------------
    // Brass
    //
    // 70% Cu + 30% Zn by mass.
    // --------------------------------------------------------

    if (!G4Material::GetMaterial("Brass", false))
    {
        auto* Cu = nist->FindOrBuildMaterial("G4_Cu");
        auto* Zn = nist->FindOrBuildMaterial("G4_Zn");

        auto* brass =
            new G4Material(
                "Brass",
                8.53 * g / cm3,
                2
            );

        brass->AddMaterial(Cu, 0.70);
        brass->AddMaterial(Zn, 0.30);
    }

    // --------------------------------------------------------
    // Plastic scintillator
    // --------------------------------------------------------

    if (!G4Material::GetMaterial("PlasticScintillator", false))
    {
        auto* H = nist->FindOrBuildElement("H");
        auto* C = nist->FindOrBuildElement("C");

        auto* scintillator =
            new G4Material(
                "PlasticScintillator",
                1.032 * g / cm3,
                2
            );

        scintillator->AddElement(C, 0.915);
        scintillator->AddElement(H, 0.085);
    }

    // --------------------------------------------------------
    // CO2
    // --------------------------------------------------------

    if (!G4Material::GetMaterial("CO2", false))
    {
        auto* C = nist->FindOrBuildElement("C");
        auto* O = nist->FindOrBuildElement("O");

        auto* co2 =
            new G4Material(
                "CO2",
                1.842e-3 * g / cm3,
                2
            );

        co2->AddElement(C, 1);
        co2->AddElement(O, 2);
    }

    // --------------------------------------------------------
    // Ar/CO2 gas mixture
    //
    // Existing simulation composition retained.
    // --------------------------------------------------------

    if (!G4Material::GetMaterial("ArCO2", false))
    {
        auto* Ar = nist->FindOrBuildElement("Ar");
        auto* C  = nist->FindOrBuildElement("C");
        auto* O  = nist->FindOrBuildElement("O");

        auto* arCO2 =
            new G4Material(
                "ArCO2",
                1.822e-3 * g / cm3,
                3
            );

        arCO2->AddElement(Ar, 0.7841);
        arCO2->AddElement(C,  0.0589);
        arCO2->AddElement(O,  0.1570);
    }
}


// ============================================================
// Construct detector
// ============================================================

G4VPhysicalVolume*
DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    auto* air =
        nist->FindOrBuildMaterial("G4_AIR");

    auto* leadGlass =
        G4Material::GetMaterial("LeadGlass");

    auto* scintillator =
        G4Material::GetMaterial("PlasticScintillator");

    auto* co2 =
        G4Material::GetMaterial("CO2");

    auto* arCO2 =
        G4Material::GetMaterial("ArCO2");

    // --------------------------------------------------------
    // World
    // --------------------------------------------------------

    const G4double worldSizeXY = 60.0 * cm;
    const G4double worldSizeZ  = 460.0 * cm;

    auto* worldSolid =
        new G4Box(
            "World",
            worldSizeXY / 2.0,
            worldSizeXY / 2.0,
            worldSizeZ / 2.0
        );

    auto* worldLogical =
        new G4LogicalVolume(
            worldSolid,
            air,
            "World"
        );

    auto* worldPhysical =
        new G4PVPlacement(
            nullptr,
            G4ThreeVector(),
            worldLogical,
            "World",
            nullptr,
            false,
            0,
            true
        );

    // --------------------------------------------------------
    // Cherenkov detector windows / gas
    //
    // Retained from the existing geometry.
    // --------------------------------------------------------

    auto* cherenkovWindowSolid =
        new G4Box(
            "CherenkovWindow",
            10.0 * cm,
            10.0 * cm,
            0.5 * mm
        );

    auto* cherenkovWindowLogical =
        new G4LogicalVolume(
            cherenkovWindowSolid,
            air,
            "CherenkovWindow"
        );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, -100.0 * cm),
        cherenkovWindowLogical,
        "CherenkovWindow",
        worldLogical,
        false,
        0,
        true
    );

    auto* cherenkovGasSolid =
        new G4Box(
            "CherenkovGas",
            10.0 * cm,
            10.0 * cm,
            50.0 * cm
        );

    auto* cherenkovGasLogical =
        new G4LogicalVolume(
            cherenkovGasSolid,
            arCO2,
            "CherenkovGas"
        );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, -50.0 * cm),
        cherenkovGasLogical,
        "CherenkovGas",
        worldLogical,
        false,
        0,
        true
    );

    // --------------------------------------------------------
    // Drift chambers
    // --------------------------------------------------------

    auto* dwcSolid =
        new G4Box(
            "DWC",
            10.0 * cm,
            10.0 * cm,
            1.0 * mm
        );

    auto* dwcLogical =
        new G4LogicalVolume(
            dwcSolid,
            arCO2,
            "DWC"
        );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, -85.0 * cm),
        dwcLogical,
        "DWC1",
        worldLogical,
        false,
        0,
        true
    );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, -70.0 * cm),
        dwcLogical,
        "DWC2",
        worldLogical,
        false,
        1,
        true
    );

    // --------------------------------------------------------
    // Scintillator 1
    // --------------------------------------------------------

    auto* scintillatorSolid =
        new G4Box(
            "Scintillator",
            10.0 * cm,
            10.0 * cm,
            2.5 * mm
        );

    auto* scintillatorLogical =
        new G4LogicalVolume(
            scintillatorSolid,
            scintillator,
            "Scintillator"
        );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, -40.0 * cm),
        scintillatorLogical,
        "Scintillator1",
        worldLogical,
        false,
        0,
        true
    );

    // --------------------------------------------------------
    // Absorber
    //
    // Downstream face remains fixed at -6 cm.
    // Therefore its center moves when thickness changes.
    // --------------------------------------------------------

    const G4double absorberCenterZ =
        kAbsorberDownstreamZ -
        fAbsorberThickness / 2.0;

    auto* absorberSolid =
        new G4Box(
            "Absorber",
            kAbsorberSizeXY / 2.0,
            kAbsorberSizeXY / 2.0,
            fAbsorberThickness / 2.0
        );

    fAbsorberLogical =
        new G4LogicalVolume(
            absorberSolid,
            fAbsorberMaterial,
            "Absorber"
        );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, absorberCenterZ),
        fAbsorberLogical,
        "Absorber",
        worldLogical,
        false,
        0,
        true
    );

    // --------------------------------------------------------
    // Scintillator 2
    // --------------------------------------------------------

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, -5.0 * cm),
        scintillatorLogical,
        "Scintillator2",
        worldLogical,
        false,
        1,
        true
    );

    // --------------------------------------------------------
    // Lead-glass calorimeter
    //
    // 4 x 4 array
    // Each block = 10 x 10 x 37 cm
    // --------------------------------------------------------

    auto* blockSolid =
        new G4Box(
            "LeadGlassBlock",
            kBlockSizeXY / 2.0,
            kBlockSizeXY / 2.0,
            kBlockSizeZ / 2.0
        );

    const G4double firstX =
        -(kNCols - 1) * kBlockSizeXY / 2.0;

    const G4double firstY =
        -(kNRows - 1) * kBlockSizeXY / 2.0;

    for (G4int row = 0; row < kNRows; ++row)
    {
        for (G4int col = 0; col < kNCols; ++col)
        {
            const G4double x =
                firstX + col * kBlockSizeXY;

            const G4double y =
                firstY + row * kBlockSizeXY;

            const G4int copyNumber =
                row * kNCols + col;

            // Unique logical-volume names make the sensitive
            // detector assignment robust.
            const G4String logicalName =
                "Block_" +
                std::to_string(row) +
                "_" +
                std::to_string(col);

            auto* blockLogical =
                new G4LogicalVolume(
                    blockSolid,
                    leadGlass,
                    logicalName
                );

            new G4PVPlacement(
                nullptr,
                G4ThreeVector(x, y, kBlockSizeZ / 2.0),
                blockLogical,
                logicalName,
                worldLogical,
                false,
                copyNumber,
                true
            );
        }
    }

    fScoringVolume =
        G4LogicalVolumeStore::GetInstance()
            ->GetVolume("Block_1_1");

    // --------------------------------------------------------
    // Visualisation
    // --------------------------------------------------------

    worldLogical->SetVisAttributes(
        G4VisAttributes::GetInvisible()
    );

    auto* absorberVis =
        new G4VisAttributes(
            G4Colour(0.7, 0.7, 0.7)
        );

    absorberVis->SetForceSolid(true);

    fAbsorberLogical->SetVisAttributes(
        absorberVis
    );

    auto* blockVis =
        new G4VisAttributes(
            G4Colour(0.2, 0.8, 1.0)
        );

    blockVis->SetForceSolid(true);

    for (auto* logicalVolume :
         *G4LogicalVolumeStore::GetInstance())
    {
        if (logicalVolume->GetName().find("Block_") == 0)
        {
            logicalVolume->SetVisAttributes(blockVis);
        }
    }

    PrintParameters();

    return worldPhysical;
}


// ============================================================
// Sensitive detector
// ============================================================

void DetectorConstruction::ConstructSDandField()
{
    auto* sdManager =
        G4SDManager::GetSDMpointer();

    auto* calorimeterSD =
        static_cast<CalorimeterSD*>(
            sdManager->FindSensitiveDetector(
                "CalorimeterSD",
                false
            )
        );

    if (!calorimeterSD)
    {
        calorimeterSD =
            new CalorimeterSD(
                "CalorimeterSD",
                "CaloHits",
                kNCols * kNRows
            );

        sdManager->AddNewDetector(
            calorimeterSD
        );
    }

    auto* logicalVolumeStore =
        G4LogicalVolumeStore::GetInstance();

    for (auto* logicalVolume :
         *logicalVolumeStore)
    {
        const G4String& name =
            logicalVolume->GetName();

        if (name.find("Block_") == 0)
        {
            logicalVolume->SetSensitiveDetector(
                calorimeterSD
            );
        }
    }
}


// ============================================================
// Set absorber material
// ============================================================

void DetectorConstruction::SetAbsorberMaterial(
    const G4String& materialName
)
{
    auto* material =
        G4Material::GetMaterial(
            materialName,
            false
        );

    if (!material)
    {
        G4Exception(
            "DetectorConstruction::SetAbsorberMaterial",
            "DetectorMaterial001",
            JustWarning,
            ("Material '" + materialName +
             "' was not found.").c_str()
        );

        return;
    }

    fAbsorberMaterial = material;

    if (fAbsorberLogical)
    {
        fAbsorberLogical->SetMaterial(
            fAbsorberMaterial
        );
    }

    G4RunManager::GetRunManager()
        ->PhysicsHasBeenModified();

    PrintParameters();
}


// ============================================================
// Set absorber thickness
// ============================================================

void DetectorConstruction::SetAbsorberThickness(
    G4double thickness
)
{
    if (thickness <= 0.0)
    {
        G4Exception(
            "DetectorConstruction::SetAbsorberThickness",
            "DetectorThickness001",
            JustWarning,
            "Absorber thickness must be greater than zero."
        );

        return;
    }

    fAbsorberThickness = thickness;

    UpdateGeometry();
}


// ============================================================
// Update geometry
// ============================================================

void DetectorConstruction::UpdateGeometry()
{
    auto* runManager =
        G4RunManager::GetRunManager();

    runManager->ReinitializeGeometry();

    PrintParameters();
}


// ============================================================
// Accessors
// ============================================================

G4Material*
DetectorConstruction::GetAbsorberMaterial() const
{
    return fAbsorberMaterial;
}


G4double
DetectorConstruction::GetAbsorberThickness() const
{
    return fAbsorberThickness;
}


// ============================================================
// Print detector parameters
// ============================================================

void DetectorConstruction::PrintParameters() const
{
    if (!fAbsorberMaterial)
    {
        return;
    }

    G4cout
        << "\n"
        << "========================================\n"
        << " Detector configuration\n"
        << "========================================\n"
        << " Absorber material : "
        << fAbsorberMaterial->GetName()
        << "\n"
        << " Absorber thickness: "
        << fAbsorberThickness / mm
        << " mm\n"
        << " Radiation length  : "
        << fAbsorberMaterial->GetRadlen() / mm
        << " mm\n"
        << " Thickness / X0    : "
        << fAbsorberThickness /
           fAbsorberMaterial->GetRadlen()
        << "\n"
        << "========================================\n"
        << G4endl;
}

