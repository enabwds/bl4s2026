#ifndef ACTION_INITIALISATION_HH
#define ACTION_INITIALISATION_HH

#include "G4VUserActionInitialization.hh"

class ActionInitialisation : public G4VUserActionInitialization
{
public:
    ActionInitialisation() = default;
    ~ActionInitialisation() override = default;

    void BuildForMaster() const override;
    void Build() const override;
};

#endif