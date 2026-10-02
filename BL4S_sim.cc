// ============================================================
//  BL4S_sim.cc
//  Main entry point for the GEANT4 electromagnetic-shower simulator
// ============================================================

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4ios.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialisation.hh"

int main(int argc, char** argv)
{
    // --------------------------------------------------------
    // 1. Create the GEANT4 run manager
    // --------------------------------------------------------
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

    // --------------------------------------------------------
    // 2. Register detector, physics, and user actions
    // --------------------------------------------------------
    runManager->SetUserInitialization(
        new DetectorConstruction()
    );

    G4PhysListFactory factory;
    G4VModularPhysicsList* physicsList =
        factory.GetReferencePhysList("FTFP_BERT_EMZ");

    G4cout << "Physics list: FTFP_BERT_EMZ" << G4endl;

    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(
        new ActionInitialisation()
    );

    // --------------------------------------------------------
    // 3. Initialise visualisation infrastructure
    //
    // We initialise the visualisation manager here so that
    // visualisation commands are available if requested.
    //
    // We DO NOT automatically open a graphics window.
    // --------------------------------------------------------
    auto* visManager = new G4VisExecutive();
    visManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();

    // --------------------------------------------------------
    // 4. Batch or interactive mode
    // --------------------------------------------------------
    if (argc > 1)
    {
        // Batch mode:
        //   ./BL4S_sim macros/test_beam.mac
        //
        // The supplied macro controls initialization and runs.
        G4String command = "/control/execute ";
        G4String fileName = argv[1];

        UI->ApplyCommand(command + fileName);
    }
    else
    {
        // Interactive mode:
        // Start the GEANT4 command interface directly.
        //
        // No visualization macro is automatically executed.
        // You should get the Idle> prompt.
        auto* ui = new G4UIExecutive(argc, argv);

        ui->SessionStart();

        delete ui;
    }

    // --------------------------------------------------------
    // 5. Clean up
    // --------------------------------------------------------
    delete visManager;
    delete runManager;

    return 0;
}