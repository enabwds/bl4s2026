#pragma once
// ============================================================
//  ActionInitialisation.hh
//  Registers all user action classes with the run manager.
// ============================================================

#include "G4VUserActionInitialization.hh"

class ActionInitialisation : public G4VUserActionInitialization
{
public:
    ActionInitialisation() = default;
    ~ActionInitialisation() override = default;

    void BuildForMaster() const override;  // MT: master thread
    void Build()          const override;  // MT: worker threads / serial
};
