// ============================================================
//  DetectorConstruction.cc
// ============================================================
//  GEOMETRY OVERVIEW (beam travels along +Z):
//
//  [Cherenkov 1] [Cherenkov 2] [DWC 1] [DWC 2]
//       z=-300        z=-250    z=-200   z=-150
//
//  [Scint 1]  [Absorber slab]  [Scint 2]  [4×4 PbGlass array]
//    z=-100       z=-50          z=-10        z=0 → z=37cm
//
//  Real BL4S block dimensions: 10×10×37 cm each, 4×4 = 40×40 cm face
//  Beam spot: 2 cm diameter circular cross-section
//  Beam divergence: ~1 mrad
// ============================================================

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "CalorimeterSD.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ios.hh"

// ── Constructor ──────────────────────────────────────────────────────────────
DetectorConstruction::DetectorConstruction()
{
    fMessenger = new DetectorMessenger(this);
    DefineMaterials();
    fAbsorberMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
}

// ── Material definitions ─────────────────────────────────────────────────────
void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

    // ── Standard absorber materials ──────────────────────────────────
    //  NIST name   |  Z  |  ρ (g/cm³) | X0 (cm) | Ec (MeV) | RM (cm)
    //  G4_Al       |  13 |   2.70     |  8.897  |  42.3    |  4.4
    //  G4_Fe       |  26 |   7.87     |  1.757  |  21.2    |  1.72
    //  G4_Cu       |  29 |   8.96     |  1.436  |  19.0    |  1.57
    //  G4_Pb       |  82 | 11.35      |  0.5612 |   6.71   |  1.60
    nist->FindOrBuildMaterial("G4_Al");
    nist->FindOrBuildMaterial("G4_Fe");
    nist->FindOrBuildMaterial("G4_Cu");
    nist->FindOrBuildMaterial("G4_Pb");
    nist->FindOrBuildMaterial("G4_AIR");
    // ── Lead glass (BL4S calorimeter blocks) ─────────────────────────
    // Composition: ~65% SiO2, ~35% PbO by mass (SF5 type lead glass)
    // Density: 3.86 g/cm³, X0 ≈ 2.37 cm, ~18 X0 in 37 cm depth
    G4Element* Si = nist->FindOrBuildElement("Si");
    G4Element* Pb = nist->FindOrBuildElement("Pb");
    G4Element* Ox = nist->FindOrBuildElement("O");
    auto* leadGlass = new G4Material("LeadGlass", 3.86*g/cm3, 3);
    leadGlass->AddElement(Pb, 0.214);   // mass fraction
    leadGlass->AddElement(Si, 0.192);
    leadGlass->AddElement(Ox, 0.594);

    // ── Brass (blind test material: 70% Cu, 30% Zn by mass) ─────────
    // Effective Z ≈ 27.4, effective X0 ≈ 1.54 cm
    // Sits between Fe and Cu — a fair but genuine inverse-problem challenge
    G4Element* Cu = nist->FindOrBuildElement("Cu");
    G4Element* Zn = nist->FindOrBuildElement("Zn");
    G4double brassDensity = 8.53 * g/cm3;
    auto* brass = new G4Material("Brass", brassDensity, 2);
    brass->AddElement(Cu, 0.70);
    brass->AddElement(Zn, 0.30);

    // ── Tungsten (harder blind test — extrapolates beyond training range)
    // Z=74, X0=0.35 cm
    nist->FindOrBuildMaterial("G4_W");

    // ── Scintillator material (plastic, for trigger scintillators) ───
    nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    // ── Cherenkov gas (CO2 at ~1 atm, pressure-tunable in real life) ─
    nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");

    // ── DWC gas (Ar/CO2 80:20 by volume) ────────────────────────────
    G4Element* Ar = nist->FindOrBuildElement("Ar");
    G4Element* C  = nist->FindOrBuildElement("C");
    G4Element* O  = nist->FindOrBuildElement("O");
    // Approximate ArCO2 as a mixture
    G4double arco2Density = 1.72e-3 * g/cm3;  // at STP
    auto* arco2 = new G4Material("ArCO2_80_20", arco2Density, 3);
    arco2->AddElement(Ar, 0.783);
    arco2->AddElement(C,  0.059);
    arco2->AddElement(O,  0.158);

    (void)brass; (void)arco2;  // suppress unused warnings
}

// ── Geometry construction ────────────────────────────────────────────────────
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Close geometry first if open, then clean all stores.
    // This must happen before any new volumes are created, and prevents
    // stale LVs from accumulating in the store across geometry rebuilds.
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4SolidStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* air     = nist->FindOrBuildMaterial("G4_AIR");

    // ── WORLD ────────────────────────────────────────────────────────
    // Large enough to hold full beamline: 400 cm long, 60×60 cm face
    auto* worldSolid = new G4Box("World", 30.*cm, 30.*cm, 200.*cm);
    auto* worldLV    = new G4LogicalVolume(worldSolid, air, "World");
    auto* worldPV    = new G4PVPlacement(nullptr, G4ThreeVector(),
                                          worldLV, "World",
                                          nullptr, false, 0, true);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ── CHERENKOV DETECTORS (1 & 2) ──────────────────────────────────
    // Each detector: 1mm Al window + 50cm gas volume, no overlaps
    // Cherenkov 1: window at z=-190, gas z=-189 to z=-139
    // Cherenkov 2: window at z=-130, gas z=-129 to z=-79
    G4Material* co2      = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    G4Material* alMat    = nist->FindOrBuildMaterial("G4_Al");
    auto* chkvWindowSolid = new G4Box("ChkvWindow", 5.*cm, 5.*cm, 0.05*cm);
    auto* chkvGasSolid    = new G4Box("ChkvGas",    5.*cm, 5.*cm, 25.*cm);
    auto* chkvWindowLV1   = new G4LogicalVolume(chkvWindowSolid, alMat, "ChkvWindow1");
    auto* chkvGasLV1      = new G4LogicalVolume(chkvGasSolid,    co2,   "ChkvGas1");
    auto* chkvWindowLV2   = new G4LogicalVolume(chkvWindowSolid, alMat, "ChkvWindow2");
    auto* chkvGasLV2      = new G4LogicalVolume(chkvGasSolid,    co2,   "ChkvGas2");

    // Window at front, gas immediately behind (no gap, no overlap)
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-189.95*cm), chkvWindowLV1, "ChkvWindow1", worldLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-164.00*cm), chkvGasLV1,    "ChkvGas1",    worldLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-138.95*cm), chkvWindowLV2, "ChkvWindow2", worldLV, false, 1, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-113.00*cm), chkvGasLV2,    "ChkvGas2",    worldLV, false, 1, true);

    auto* chkvVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.0, 0.3));
    chkvVis->SetForceSolid(true);
    chkvWindowLV1->SetVisAttributes(chkvVis);
    chkvWindowLV2->SetVisAttributes(chkvVis);
    auto* chkvGasVis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.1));
    chkvGasLV1->SetVisAttributes(chkvGasVis);
    chkvGasLV2->SetVisAttributes(chkvGasVis);

    // ── DELAY WIRE CHAMBERS (DWC 1 & 2) ─────────────────────────────
    // 10×10 cm, 1 cm thick Ar/CO2 gas, no overlap with Cherenkov
    G4Material* arco2 = G4Material::GetMaterial("ArCO2_80_20");
    auto* dwcSolid    = new G4Box("DWC", 5.*cm, 5.*cm, 0.5*cm);
    auto* dwcLV1      = new G4LogicalVolume(dwcSolid, arco2, "DWC1");
    auto* dwcLV2      = new G4LogicalVolume(dwcSolid, arco2, "DWC2");

    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-85.*cm), dwcLV1, "DWC1", worldLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-70.*cm), dwcLV2, "DWC2", worldLV, false, 1, true);

    auto* dwcVis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.5, 0.5));
    dwcVis->SetForceSolid(true);
    dwcLV1->SetVisAttributes(dwcVis);
    dwcLV2->SetVisAttributes(dwcVis);

    // ── SCINTILLATOR TRIGGER 1 (upstream of absorber) ────────────────
    // 10×20 cm plastic scintillator, 5 mm thick
    G4Material* scint = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    auto* scint1Solid = new G4Box("Scint1", 5.*cm, 10.*cm, 0.25*cm);
    auto* scint1LV    = new G4LogicalVolume(scint1Solid, scint, "Scint1");
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-40.*cm),
                      scint1LV, "Scint1", worldLV, false, 0, true);

    auto* scintVis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.0, 0.7));
    scintVis->SetForceSolid(true);
    scint1LV->SetVisAttributes(scintVis);

    // ── ABSORBER SLAB ────────────────────────────────────────────────
    // Variable element — material and thickness set via macro commands.
    // Must be wide enough to catch full beam halo (beam spot = 2 cm diam).
    // Downstream face anchored at kAbsorberDownstreamZ so the air gap to
    // the calorimeter stays constant regardless of absorber thickness.
    G4double absorberHalfZ   = 0.5 * fAbsorberThickness * mm;
    G4double absorberCentreZ = kAbsorberDownstreamZ*cm - absorberHalfZ;
    auto* absorberSolid = new G4Box("Absorber",
        10.*cm, 10.*cm, absorberHalfZ);
    auto* absorberLV = new G4LogicalVolume(absorberSolid,
                                            fAbsorberMaterial, "Absorber");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, absorberCentreZ),
                      absorberLV, "Absorber", worldLV, false, 0, true);

    auto* absorberVis = new G4VisAttributes(G4Colour(0.8, 0.2, 0.2, 0.8));
    absorberVis->SetForceSolid(true);
    absorberLV->SetVisAttributes(absorberVis);

    // ── SCINTILLATOR TRIGGER 2 (downstream of absorber) ──────────────
    auto* scint2Solid = new G4Box("Scint2", 5.*cm, 10.*cm, 0.25*cm);
    auto* scint2LV    = new G4LogicalVolume(scint2Solid, scint, "Scint2");
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-5.*cm),
                      scint2LV, "Scint2", worldLV, false, 1, true);
    scint2LV->SetVisAttributes(scintVis);

    // ── CALORIMETER ARRAY (4×4 lead-glass blocks) ────────────────────
    // Real BL4S specs: 10×10×37 cm blocks, 4×4 array = 40×40 cm face
    // ~18 radiation lengths in lead glass — fully contains GeV EM showers
    // Energy resolution: σ_E/E = 0.02% ⊕ 6.3%/√E
    G4Material* pbGlass   = G4Material::GetMaterial("LeadGlass");
    G4double caloFrontZ   = 0.*cm;   // front face at z=0
    G4double caloHalfZ    = 0.5 * kBlockSizeZ * cm;
    G4double caloMidZ     = caloFrontZ + caloHalfZ;
    G4double arrayHalfXY  = 0.5 * kNcols * kBlockSizeXY * cm;

    auto* caloSolid = new G4Box("Calorimeter",
        arrayHalfXY, arrayHalfXY, caloHalfZ);
    auto* caloLV = new G4LogicalVolume(caloSolid, air, "Calorimeter");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, caloMidZ),
                      caloLV, "Calorimeter", worldLV, false, 0, true);

    auto* blockVis = new G4VisAttributes(G4Colour(0.2, 0.5, 1.0, 0.6));
    blockVis->SetForceSolid(true);

    // Give each block its own uniquely named logical volume.
    // This avoids the "More than one logical volume named Block" crash
    // that occurs when GEANT4 MT tries to reassign the SD after geometry
    // reinitialisation — each LV name must be unique for SetSensitiveDetector.
    auto* blockSolid = new G4Box("BlockSolid",
        0.5*kBlockSizeXY*cm, 0.5*kBlockSizeXY*cm, caloHalfZ);

    int copyNo = 0;
    G4LogicalVolume* firstBlockLV = nullptr;
    for (int row = 0; row < kNrows; ++row) {
        for (int col = 0; col < kNcols; ++col) {
            G4String lvName = "Block_" + std::to_string(row) + "_" + std::to_string(col);
            auto* blockLV = new G4LogicalVolume(blockSolid, pbGlass, lvName);
            blockLV->SetVisAttributes(blockVis);
            G4double xPos = (-1.5 + col) * kBlockSizeXY*cm;
            G4double yPos = (-1.5 + row) * kBlockSizeXY*cm;
            new G4PVPlacement(nullptr,
                              G4ThreeVector(xPos, yPos, 0),
                              blockLV, lvName,
                              caloLV, false, copyNo++, true);
            if (!firstBlockLV) firstBlockLV = blockLV;
        }
    }

    fScoringVolume = firstBlockLV;
    PrintParameters();
    return worldPV;
}

// ── Sensitive Detectors ──────────────────────────────────────────────────────
void DetectorConstruction::ConstructSDandField()
{
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    // SD is shared across threads — create once on first call, reuse thereafter.
    auto* caloSD = dynamic_cast<CalorimeterSD*>(
        sdManager->FindSensitiveDetector("CalorimeterSD", false));
    if (!caloSD) {
        caloSD = new CalorimeterSD("CalorimeterSD", "HitsCollection",
                                   kNrows * kNcols);
        sdManager->AddNewDetector(caloSD);
    }

    // Assign SD by iterating the LV store directly, matching on name prefix.
    //
    // Why not SetSensitiveDetector(name, sd)?
    //   That helper does a name lookup in G4LogicalVolumeStore. In MT mode,
    //   ConstructSDandField() runs on each worker thread, but Construct() only
    //   runs on the master. Workers clone the master geometry, so the LVs exist
    //   in the store, but any name-suffix scheme tied to fGeneration is fragile
    //   because workers may see different store state.
    //
    // Iterating the store and matching on the "Block_" prefix is robust:
    //   it finds whatever LVs were actually created, regardless of suffix,
    //   and assigns the SD directly to the pointer — no name lookup needed.
    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    int assigned = 0;
    for (auto* lv : *lvStore) {
        const G4String& name = lv->GetName();
        if (name.substr(0, 6) == "Block_") {
            lv->SetSensitiveDetector(caloSD);
            ++assigned;
        }
    }
    if (assigned == 0) {
        G4cerr << "ConstructSDandField: WARNING — no Block_ LVs found in store!\n";
    }
}

// ── Setters ──────────────────────────────────────────────────────────────────
void DetectorConstruction::SetAbsorberMaterial(const G4String& name)
{
    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(name);
    if (!mat) mat = G4Material::GetMaterial(name);
    if (!mat) { G4cerr << "Material " << name << " not found!\n"; return; }
    fAbsorberMaterial = mat;

    // Update the existing logical volume directly — no geometry rebuild needed.
    // This avoids the MT duplicate-LV crash caused by ReinitializeGeometry().
    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    auto* absorberLV = lvStore->GetVolume("Absorber", false);
    if (absorberLV) absorberLV->SetMaterial(mat);

    // Notify GEANT4 that material changed (needed for physics tables)
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetAbsorberThickness(G4double t)
{
    fAbsorberThickness = t / mm;
    // Full geometry rebuild so both the solid size AND the centre position
    // (which depends on thickness via kAbsorberDownstreamZ - halfZ) update
    // together. In-place SetZHalfLength() only resized the solid but left
    // the placement centre fixed at z=-20cm, causing the downstream face
    // to drift by up to ~4cm across the thickness scan.
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    PrintParameters();
}

// ── Diagnostics ──────────────────────────────────────────────────────────────
void DetectorConstruction::PrintParameters() const
{
    G4double X0 = fAbsorberMaterial->GetRadlen();
    G4cout << "\n──── Detector Parameters ──────────────────────────────────────\n"
           << "  Absorber material  : " << fAbsorberMaterial->GetName() << "\n"
           << "  Absorber thickness : " << fAbsorberThickness << " mm"
           << "  (" << (fAbsorberThickness*mm)/X0 << " X0)\n"
           << "  Radiation length   : " << X0/cm << " cm\n"
           << "  Calorimeter        : " << kNrows << "×" << kNcols
           << " blocks, each " << kBlockSizeXY << "×" << kBlockSizeXY
           << "×" << kBlockSizeZ << " cm\n"
           << "  Calo face area     : "
           << kNcols*kBlockSizeXY << "×" << kNrows*kBlockSizeXY << " cm\n"
           << "───────────────────────────────────────────────────────────────\n";
}
