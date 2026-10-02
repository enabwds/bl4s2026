#pragma once

// ============================================================
//  CalorimeterSD.hh
//
//  Sensitive detector for the 4x4 lead-glass calorimeter.
//
//  Responsibility:
//    - Create one hit for each calorimeter block.
//    - Accumulate ALL positive energy deposition in each block.
//    - Make the raw per-block energy available to EventAction.
//
//  This class deliberately does NOT:
//    - apply energy thresholds
//    - reject events
//    - calculate shower features
//    - perform noise filtering
//    - perform ML preprocessing
//
//  Those operations belong to the analysis layer.
// ============================================================

#include "G4VHit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"
#include "globals.hh"

// ============================================================
//  Calorimeter hit
// ============================================================

class CaloHit : public G4VHit
{
public:
    CaloHit()
        : fEdep(0.0),
          fBlockID(-1)
    {}

    ~CaloHit() override = default;

    // Accumulate energy deposited in this block.
    void AddEdep(G4double energy)
    {
        fEdep += energy;
    }

    // Return total energy deposited in this block
    // during the current event.
    G4double GetEdep() const
    {
        return fEdep;
    }

    // Block index: 0 ... 15.
    G4int GetBlockID() const
    {
        return fBlockID;
    }

    void SetBlockID(G4int blockID)
    {
        fBlockID = blockID;
    }

private:
    G4double fEdep;
    G4int    fBlockID;
};


// ============================================================
//  Hit collection
// ============================================================

using CaloHitsCollection = G4THitsCollection<CaloHit>;


// ============================================================
//  Sensitive detector
// ============================================================

class CalorimeterSD : public G4VSensitiveDetector
{
public:

    CalorimeterSD(const G4String& name,
                  const G4String& hitsCollectionName,
                  G4int nCells);

    ~CalorimeterSD() override = default;

    // Called once at the beginning of every event.
    void Initialize(G4HCofThisEvent* hce) override;

    // Called for every GEANT4 step inside a sensitive
    // calorimeter block.
    G4bool ProcessHits(G4Step* step,
                       G4TouchableHistory* history) override;

    // Called at the end of every event.
    void EndOfEvent(G4HCofThisEvent* hce) override;

private:

    // Current event's hit collection.
    CaloHitsCollection* fHitsCollection = nullptr;

    // GEANT4 collection ID.
    G4int fHCID = -1;

    // Number of calorimeter cells.
    G4int fNCells = 0;
};
