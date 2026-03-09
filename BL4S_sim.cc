// ============================================================
//  BL4S_sim.cc  –  Main entry point
//  Electromagnetic Shower / Inverse Material Characterisation
//  Beams for Schools 2026 experiment
// ============================================================
//
//  HOW TO BUILD (once GEANT4 is installed & sourced):
//    mkdir build && cd build
//    cmake ..
//    make -j$(nproc)
//
//  HOW TO RUN:
//    ./BL4S_sim macros/run_iron.mac       # batch mode
//    ./BL4S_sim                           # interactive (opens Qt/OGL window)
//
// ============================================================

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4PhysListFactory.hh"  // FTFP_BERT_EMZ — precision EM shower physics

#include "DetectorConstruction.hh"
#include "ActionInitialisation.hh"

int main(int argc, char** argv)
{
    // ── 1. Create the run manager ──────────────────────────────────────
    // G4RunManagerFactory chooses Serial / MT automatically.
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    // ── 2. Register user-defined classes ──────────────────────────────
    runManager->SetUserInitialization(new DetectorConstruction());

    // Physics list: FTFP_BERT_EMZ
    // EMZ activates G4EmStandardPhysics_option4 — the highest-accuracy EM
    // models in GEANT4, specifically validated for GeV-scale EM calorimetry.
    // This correctly models bremsstrahlung and pair production cross-sections,
    // which directly determine CoreFraction, ShowerWidth, and attenuation curves.
    // QBBC (previous list) is a general hadronic list not optimised for this.
    G4PhysListFactory factory;
    G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("FTFP_BERT_EMZ");
    G4cout << "Physics list: FTFP_BERT_EMZ\n";
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new ActionInitialisation());

    // ── 3. UI / Visualisation ─────────────────────────────────────────
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();

    if (argc > 1) {
        // Batch mode: execute the macro file given on the command line
        G4String command  = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command + fileName);
    } else {
        // Interactive mode: open a UI session
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        UI->ApplyCommand("/control/execute macros/vis.mac");
        ui->SessionStart();
        delete ui;
    }

    // ── 4. Clean up ───────────────────────────────────────────────────
    delete visManager;
    delete runManager;
    return 0;
}
