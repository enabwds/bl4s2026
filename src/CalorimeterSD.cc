// ============================================================
//  CalorimeterSD.cc
// ============================================================

#include "CalorimeterSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"


// ============================================================
//  Constructor
// ============================================================

CalorimeterSD::CalorimeterSD(const G4String& name,
                             const G4String& hitsCollectionName,
                             G4int nCells)
    : G4VSensitiveDetector(name),
      fNCells(nCells)
{
    collectionName.insert(hitsCollectionName);
}


// ============================================================
//  Initialize
//
//  Create a fresh hit collection for every event.
//
//  We pre-create one hit for every calorimeter block so that
//  EventAction can safely access:
//
//      hit[0] ... hit[15]
//
//  even when some blocks receive zero energy.
// ============================================================

void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new CaloHitsCollection(
        SensitiveDetectorName,
        collectionName[0]
    );

    // Obtain the collection ID once and reuse it.
    if (fHCID < 0)
    {
        fHCID = G4SDManager::GetSDMpointer()
                    ->GetCollectionID(collectionName[0]);
    }

    hce->AddHitsCollection(fHCID, fHitsCollection);

    // One hit object per calorimeter block.
    for (G4int blockID = 0; blockID < fNCells; ++blockID)
    {
        auto* hit = new CaloHit();

        hit->SetBlockID(blockID);

        fHitsCollection->insert(hit);
    }
}


// ============================================================
//  ProcessHits
//
//  Called for every GEANT4 step inside a sensitive block.
//
//  IMPORTANT:
//  There is intentionally NO per-step energy threshold.
//
//  Tiny deposits must be accumulated rather than discarded.
//  Any detector threshold/noise model can be applied later in
//  the analysis layer without changing the underlying GEANT4
//  truth.
// ============================================================

G4bool CalorimeterSD::ProcessHits(G4Step* step,
                                  G4TouchableHistory*)
{
    const G4double energyDeposit =
        step->GetTotalEnergyDeposit();

    // Ignore steps that deposit no energy.
    if (energyDeposit <= 0.0)
    {
        return false;
    }

    // The block copy number identifies which of the 16
    // calorimeter blocks received this energy deposit.
    const G4int blockID =
        step->GetPreStepPoint()
            ->GetTouchable()
            ->GetReplicaNumber(0);

    // Safety check before indexing the collection.
    if (blockID < 0 || blockID >= fNCells)
    {
        return false;
    }

    // Accumulate the energy deposit.
    (*fHitsCollection)[blockID]->AddEdep(energyDeposit);

    return true;
}


// ============================================================
//  EndOfEvent
//
//  Nothing is filtered or transformed here.
//
//  EventAction is responsible for reading the completed hit
//  collection and constructing the event-level output.
// ============================================================

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
    // Intentionally empty.
}