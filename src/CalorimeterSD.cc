// ============================================================
//  CalorimeterSD.cc
// ============================================================

#include "CalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

CalorimeterSD::CalorimeterSD(const G4String& name,
                             const G4String& hitsCollectionName,
                             G4int           nCells)
: G4VSensitiveDetector(name), fNCells(nCells)
{
    collectionName.insert(hitsCollectionName);
}

// Called at the start of every event – creates fresh hit objects
void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new CaloHitsCollection(
        SensitiveDetectorName,
        collectionName[0]
    );

    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()
                    ->GetCollectionID(collectionName[0]);

    hce->AddHitsCollection(fHCID, fHitsCollection);

    // Pre-allocate one hit per block so indexing is safe
    for (int i = 0; i < fNCells; ++i) {
        auto* h = new CaloHit();
        h->SetBlockID(i);
        fHitsCollection->insert(h);
    }
}

// Called for every GEANT4 step inside a registered sensitive volume
G4bool CalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4double edep = step->GetTotalEnergyDeposit();

    // Ignore steps with no energy deposition.
    // IMPORTANT:
    // Do NOT apply a threshold here.
    // Individual small deposits must be accumulated first.
    if (edep <= 0.)
        return false;

    // copyNumber encodes which block (0–15) was hit
    G4int blockID = step->GetPreStepPoint()
                        ->GetTouchable()
                        ->GetReplicaNumber(0);

    (*fHitsCollection)[blockID]->AddEdep(edep);

    return true;
}

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
    // Optional: print a per-event summary (uncomment for debugging)
    /*
    G4cout << "  -- Calorimeter hits (>0 MeV) --\n";
    for (int i = 0; i < fNCells; ++i) {
        G4double e = (*fHitsCollection)[i]->GetEdep();
        if (e > 0.)
            G4cout << "    Block " << i << ": "
                   << e/MeV << " MeV\n";
    }
    */
}