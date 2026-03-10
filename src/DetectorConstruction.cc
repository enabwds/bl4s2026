// ============================================================
//  DetectorConstruction.cc
// ============================================================
//  GEOMETRY OVERVIEW (beam travels along +Z):
//
//  [ChkvWin1]  [ChkvGas1]       [ChkvWin2]  [ChkvGas2]
//   z=-190cm   centre=-164.90    z=-139.85   centre=-114.80
//   (flush, no gaps between window and gas volumes)
//
//  [DWC 1]  [DWC 2]  [Scint 1]  [Absorber slab]  [Scint 2]  [4×4 PbGlass array]
//   z=-85cm   z=-70cm   z=-40cm   downstream       z=-5cm      z=0 → z=+37cm
//                                 face at z=-6cm
//
//  Real BL4S block dimensions: 10×10×37 cm each, 4×4 = 40×40 cm face
//  Beam spot: 2 cm diameter circular cross-section
//  Beam divergence: ~1 mrad
//  Beam gun origin: z = -191 cm (1 cm upstream of ChkvWindow1 front face)
// ============================================================

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "CalorimeterSD.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
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
    //  G4_Pb       |  82 | 11.35      |  0.5612 |   7.43   |  1.60  (Ec: PDG 2023 electron value)
    nist->FindOrBuildMaterial("G4_Al");
    nist->FindOrBuildMaterial("G4_Fe");
    nist->FindOrBuildMaterial("G4_Cu");
    nist->FindOrBuildMaterial("G4_Pb");
    nist->FindOrBuildMaterial("G4_AIR");
    // ── Lead glass (BL4S calorimeter blocks) ─────────────────────────
    // SF5-type lead glass: ~51% PbO + ~49% SiO2 by mass.
    // Converting oxide fractions to elemental mass fractions:
    //   PbO (51%): Pb = 0.51 * (207.2/223.2) = 0.474
    //              O  = 0.51 * ( 16.0/223.2) = 0.037
    //   SiO2(49%): Si = 0.49 * ( 28.1/60.1)  = 0.229
    //              O  = 0.49 * ( 32.0/60.1)  = 0.261
    //   Total O = 0.037 + 0.261 = 0.298 (rounded so fractions sum to 1.000)
    // Density: 3.86 g/cm3.
    // Published SF5 X0 ~ 2.39-2.46 cm (density 3.86 g/cm3), giving ~15-16 X0
    // in 37 cm depth — sufficient to contain GeV-scale EM showers (need ~20 X0
    // for full containment; ~1-2% leakage expected at 4 GeV).
    // GEANT4 computes X0 internally via the Tsai formula; the simple Bragg
    // additivity rule underestimates X0 for high-Z mixtures so do not use it
    // to cross-check this value.
    G4Element* Si = nist->FindOrBuildElement("Si");
    G4Element* Pb = nist->FindOrBuildElement("Pb");
    G4Element* Ox = nist->FindOrBuildElement("O");
    auto* leadGlass = new G4Material("LeadGlass", 3.86*g/cm3, 3);
    leadGlass->AddElement(Pb, 0.474);   // mass fraction — SF5 spec
    leadGlass->AddElement(Si, 0.229);
    leadGlass->AddElement(Ox, 0.297);   // sums to 1.000

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
    // Ar/CO2 80:20 by volume at STP.
    // rho = 0.80*rho_Ar + 0.20*rho_CO2 = 0.80*1.784e-3 + 0.20*1.977e-3
    //     = 1.427e-3 + 0.395e-3 = 1.822e-3 g/cm3
    // Mass fractions from mole fractions (80% Ar, 20% CO2):
    //   M_Ar=39.948, M_CO2=44.009 -> M_mix = 0.80*39.948 + 0.20*44.009 = 40.760
    //   w_Ar = 0.80*39.948/40.760 = 0.7841
    //   w_C  = 0.20*12.011/40.760 = 0.0589
    //   w_O  = 0.20*31.998/40.760 = 0.1570  (sum = 1.0000)
    G4double arco2Density = 1.822e-3 * g/cm3;  // at STP (corrected)
    auto* arco2 = new G4Material("ArCO2_80_20", arco2Density, 3);
    arco2->AddElement(Ar, 0.7841);
    arco2->AddElement(C,  0.0589);
    arco2->AddElement(O,  0.1570);

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
    // Large enough to hold full beamline including upstream beam start
    // and full calorimeter depth: 460 cm long, 60×60 cm face
    auto* worldSolid = new G4Box("World", 30.*cm, 30.*cm, 230.*cm);
    auto* worldLV    = new G4LogicalVolume(worldSolid, air, "World");
    auto* worldPV    = new G4PVPlacement(nullptr, G4ThreeVector(),
                                          worldLV, "World",
                                          nullptr, false, 0, true);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ── CHERENKOV DETECTORS (1 & 2) ──────────────────────────────────
    // Each detector: 1 mm (0.1 cm) Al window + 50 cm gas volume.
    // Volumes are placed flush — no gaps, no overlaps.
    //
    // Layout (all values = centre position, half-extents in parentheses):
    //   ChkvWindow1: centre=-189.95 cm (halfZ=0.05) -> front=-190.00, back=-189.90
    //   ChkvGas1:    centre=-164.90 cm (halfZ=25.0)  -> front=-189.90, back=-139.90
    //   ChkvWindow2: centre=-139.85 cm (halfZ=0.05) -> front=-139.90, back=-139.80
    //   ChkvGas2:    centre=-114.80 cm (halfZ=25.0)  -> front=-139.80, back= -89.80
    G4Material* co2      = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    G4Material* alMat    = nist->FindOrBuildMaterial("G4_Al");
    auto* chkvWindowSolid = new G4Box("ChkvWindow", 5.*cm, 5.*cm, 0.05*cm);
    auto* chkvGasSolid    = new G4Box("ChkvGas",    5.*cm, 5.*cm, 25.*cm);
    auto* chkvWindowLV1   = new G4LogicalVolume(chkvWindowSolid, alMat, "ChkvWindow1");
    auto* chkvGasLV1      = new G4LogicalVolume(chkvGasSolid,    co2,   "ChkvGas1");
    auto* chkvWindowLV2   = new G4LogicalVolume(chkvWindowSolid, alMat, "ChkvWindow2");
    auto* chkvGasLV2      = new G4LogicalVolume(chkvGasSolid,    co2,   "ChkvGas2");

    // Centres calculated so each volume's back face = next volume's front face.
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-189.95*cm), chkvWindowLV1, "ChkvWindow1", worldLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-164.90*cm), chkvGasLV1,    "ChkvGas1",    worldLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-139.85*cm), chkvWindowLV2, "ChkvWindow2", worldLV, false, 1, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,-114.80*cm), chkvGasLV2,    "ChkvGas2",    worldLV, false, 1, true);

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
    // fAbsorberThickness is stored in GEANT4 internal units (mm=1).
    G4double absorberHalfZ   = 0.5 * fAbsorberThickness;
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
    // Published SF5 lead glass X0 ~2.39-2.46 cm -> ~15-16 X0 in 37 cm depth.
    // Containment is adequate for 1-2 GeV; ~1-2% longitudinal leakage at 4 GeV.
    // Energy resolution: sigma_E/E = 0.02% + 6.3%/sqrt(E)
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
           << "  Absorber thickness : " << fAbsorberThickness/mm << " mm"
           << "  (" << fAbsorberThickness/X0 << " X0)\n"
           << "  Radiation length   : " << X0/cm << " cm\n"
           << "  Calorimeter        : " << kNrows << "x" << kNcols
           << " blocks, each " << kBlockSizeXY << "x" << kBlockSizeXY
           << "x" << kBlockSizeZ << " cm\n"
           << "  Calo face area     : "
           << kNcols*kBlockSizeXY << "x" << kNrows*kBlockSizeXY << " cm\n"
           << "───────────────────────────────────────────────────────────────\n";
}
