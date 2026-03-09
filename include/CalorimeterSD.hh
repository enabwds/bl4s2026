#pragma once
// ============================================================
//  CalorimeterSD.hh
//  Sensitive Detector for the 4×4 lead-glass calorimeter.
//
//  Each block accumulates the total energy deposited by all
//  secondaries that pass through it during one event.
//  At the end of each event the hit collection is written to
//  the G4HCofThisEvent and picked up by EventAction.
// ============================================================

#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// ── Simple hit class ─────────────────────────────────────────────────────────
class CaloHit : public G4VHit
{
public:
    CaloHit() : fEdep(0.), fBlockID(-1) {}

    void AddEdep(G4double e) { fEdep += e; }
    G4double    GetEdep()    const { return fEdep; }
    G4int       GetBlockID() const { return fBlockID; }
    void        SetBlockID(G4int id) { fBlockID = id; }

private:
    G4double fEdep;
    G4int    fBlockID;
};

using CaloHitsCollection = G4THitsCollection<CaloHit>;

// ── Sensitive Detector ────────────────────────────────────────────────────────
class CalorimeterSD : public G4VSensitiveDetector
{
public:
    CalorimeterSD(const G4String& name,
                  const G4String& hitsCollectionName,
                  G4int           nCells);
    ~CalorimeterSD() override = default;

    void   Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;
    void   EndOfEvent(G4HCofThisEvent* hce) override;

private:
    CaloHitsCollection* fHitsCollection = nullptr;
    G4int               fHCID          = -1;
    G4int               fNCells;
};
