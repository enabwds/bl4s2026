#include "ActionInitialisation.hh"

#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"

// GEANT4 action classes
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

// GEANT4 basics
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"
#include "G4THitsCollection.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

// Standard C++
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <mutex>
#include <cmath>
#include <algorithm>

// ============================================================
// Configuration
// ============================================================

namespace SimulationConfig
{
    constexpr G4double kBeamMomentumSpreadFWHM = 0.15 * GeV;
    constexpr G4double kBeamAngularSigma = 0.5 * mrad;
    constexpr G4double kBeamSpotRadius = 1.0 * cm;

    constexpr G4double kBeamVertexZ = -191.0 * cm;

    constexpr G4int kNumberOfCalorimeterBlocks = 16;
}

// ============================================================
// Primary Generator
// ============================================================

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

    PrimaryGeneratorAction()
    {
        fParticleGun = new G4ParticleGun(1);

        auto* particleTable = G4ParticleTable::GetParticleTable();
        auto* electron = particleTable->FindParticle("e-");

        fParticleGun->SetParticleDefinition(electron);
        fParticleGun->SetParticleEnergy(2.0 * GeV);
        fParticleGun->SetParticlePosition(
            G4ThreeVector(0.0, 0.0, SimulationConfig::kBeamVertexZ)
        );
    }

    ~PrimaryGeneratorAction() override
    {
        delete fParticleGun;
    }

    void GeneratePrimaries(G4Event* event) override
    {
        const G4double nominalEnergy =
            fParticleGun->GetParticleEnergy();

        // ----------------------------------------------------
        // Beam momentum spread
        //
        // The supplied beam spread is interpreted as FWHM.
        // Convert FWHM -> Gaussian sigma.
        // ----------------------------------------------------

        const G4double sigmaP =
            SimulationConfig::kBeamMomentumSpreadFWHM / 2.354820045;

        G4double momentum =
            G4RandGauss::shoot(nominalEnergy, sigmaP);

        // Avoid unphysical negative/very small momentum.
        momentum = std::max(momentum, 0.1 * GeV);

        // ----------------------------------------------------
        // Beam angular divergence
        // ----------------------------------------------------

        const G4double thetaX =
            G4RandGauss::shoot(
                0.0,
                SimulationConfig::kBeamAngularSigma
            );

        const G4double thetaY =
            G4RandGauss::shoot(
                0.0,
                SimulationConfig::kBeamAngularSigma
            );

        G4ThreeVector direction(
            std::sin(thetaX),
            std::sin(thetaY),
            1.0
        );

        direction = direction.unit();

        // ----------------------------------------------------
        // Circular beam spot
        // ----------------------------------------------------

        const G4double radius =
            SimulationConfig::kBeamSpotRadius *
            std::sqrt(G4UniformRand());

        const G4double phi =
            CLHEP::twopi * G4UniformRand();

        const G4double x =
            radius * std::cos(phi);

        const G4double y =
            radius * std::sin(phi);

        // ----------------------------------------------------
        // Configure particle
        // ----------------------------------------------------

        fParticleGun->SetParticleMomentum(momentum);
        fParticleGun->SetParticleMomentumDirection(direction);

        fParticleGun->SetParticlePosition(
            G4ThreeVector(
                x,
                y,
                SimulationConfig::kBeamVertexZ
            )
        );

        fParticleGun->GeneratePrimaryVertex(event);
    }

private:

    G4ParticleGun* fParticleGun = nullptr;
};


// ============================================================
// Raw event output
// ============================================================
//
// This output is intentionally analysis-agnostic.
//
// Every generated event is written.
// No ML filtering.
// No outlier rejection.
// No ±20% beam-energy selection.
// No target-specific preprocessing.
//
// Python/analysis code decides later what constitutes a
// usable event for a particular scientific question.
// ============================================================

class CsvOutput
{
public:

    static CsvOutput& Instance()
    {
        static CsvOutput instance;
        return instance;
    }

    void Open(const std::string& filename)
    {
        std::lock_guard<std::mutex> lock(fMutex);

        fFile.open(
            filename,
            std::ios::out |
            std::ios::trunc
        );

        if (!fFile.is_open())
        {
            G4Exception(
                "CsvOutput::Open",
                "Output001",
                FatalException,
                "Could not open shower output file."
            );
        }

        WriteHeader();
    }

    void Close()
    {
        std::lock_guard<std::mutex> lock(fMutex);

        if (fFile.is_open())
        {
            fFile.flush();
            fFile.close();
        }
    }

    void WriteEvent(
    G4int runID,
    G4int eventID,
    G4double nominalBeamEnergy,
    const std::string& material,
    G4double absorberThickness,
    G4double absorberX0,
    G4double beamEnergy,
    const G4ThreeVector& beamMomentum,
    G4double vertexX,
    G4double vertexY,
    const std::vector<G4double>& blockEnergy,
    G4double totalEnergy,
    G4double coreFraction,
    G4double showerWidth
    )


    {
        std::lock_guard<std::mutex> lock(fMutex);

        if (!fFile.is_open())
        {
            return;
        }

        fFile
            << runID << ','
            << eventID << ','
            << nominalBeamEnergy / GeV << ','
            << material << ','
            << absorberThickness / mm << ','
            << absorberX0 / mm << ','
            << beamEnergy / GeV << ','
            << beamMomentum.x() / MeV << ','
            << beamMomentum.y() / MeV << ','
            << beamMomentum.z() / MeV << ','
            << vertexX / mm << ','
            << vertexY / mm;

        for (G4double energy : blockEnergy)
        {
            fFile << ',' << energy / GeV;
        }

        fFile
            << ',' << totalEnergy / GeV
            << ',' << coreFraction
            << ',' << showerWidth / cm
            << '\n';
    }

private:

    CsvOutput() = default;

    void WriteHeader()
    {
        fFile
            << "run_id,"
            << "event_id,"
            << "nominal_beam_energy_GeV,"
            << "material,"
            << "absorber_thickness_mm,"
            << "absorber_x0_mm,"
            << "beam_energy_GeV,"
            << "beam_px_MeVc,"
            << "beam_py_MeVc,"
            << "beam_pz_MeVc,"
            << "vertex_x_mm,"
            << "vertex_y_mm";

        for (G4int i = 0;
             i < SimulationConfig::kNumberOfCalorimeterBlocks;
             ++i)
        {
            fFile
                << ",E_block_"
                << i
                << "_GeV";
        }

        fFile
            << ",Etotal_GeV"
            << ",CoreFraction"
            << ",ShowerWidth_cm"
            << '\n';
    }

    std::ofstream fFile;
    std::mutex fMutex;
};


// ============================================================
// Run Action
// ============================================================

class RunAction : public G4UserRunAction
{
public:

    RunAction() = default;

    ~RunAction() override
    {
        CsvOutput::Instance().Close();
    }

    void BeginOfRunAction(const G4Run* run) override
    {
        const G4int runID = run->GetRunID();

        // Each run gets its own output file.
        const std::string filename =
            "shower_output_run_" +
            std::to_string(runID) +
            ".csv";

        CsvOutput::Instance().Open(filename);
    }

    void EndOfRunAction(const G4Run*) override
    {
        CsvOutput::Instance().Close();
    }
};


// ============================================================
// Event Action
// ============================================================

class EventAction : public G4UserEventAction
{
public:

    explicit EventAction(
        const DetectorConstruction* detector
    )
        : fDetector(detector)
    {
    }

    ~EventAction() override = default;

    void EndOfEventAction(
        const G4Event* event
    ) override
    {
        // ----------------------------------------------------
        // Retrieve calorimeter hit collection
        // ----------------------------------------------------

        auto* hce = event->GetHCofThisEvent();

        if (!hce)
        {
            return;
        }

        auto* hits =
            static_cast<G4THitsCollection<CaloHit>*>(
                hce->GetHC(
                    G4SDManager::GetSDMpointer()
                        ->GetCollectionID("CalorimeterSD/CaloHits")
                )
            );

        if (!hits)
        {
            return;
        }

        // ----------------------------------------------------
        // Read block energies
        // ----------------------------------------------------

        std::vector<G4double> blockEnergy(
            SimulationConfig::kNumberOfCalorimeterBlocks,
            0.0
        );

        G4double totalEnergy = 0.0;

        for (G4int i = 0;
             i < SimulationConfig::kNumberOfCalorimeterBlocks;
             ++i)
        {
            if ((*hits)[i])
            {
                blockEnergy[i] =
                    (*hits)[i]->GetEdep();

                totalEnergy += blockEnergy[i];
            }
        }

        // ----------------------------------------------------
        // Primary beam information
        // ----------------------------------------------------

        const auto* primaryVertex =
            event->GetPrimaryVertex();

        if (!primaryVertex)
        {
            return;
        }

        const auto* primary =
            primaryVertex->GetPrimary();

        if (!primary)
        {
            return;
        }

        const G4double beamEnergy =
            primary->GetKineticEnergy();

        const G4ThreeVector beamMomentum =
            primary->GetMomentum();

        const G4double vertexX =
            primaryVertex->GetX0();

        const G4double vertexY =
            primaryVertex->GetY0();

        // ----------------------------------------------------
        // Detector truth/configuration
        // ----------------------------------------------------

        const std::string material =
            fDetector->GetAbsorberMaterial()->GetName();

        const G4double absorberThickness =
            fDetector->GetAbsorberThickness();

        const G4double absorberX0 =
            fDetector->GetAbsorberMaterial()->GetRadlen();

        const G4double nominalBeamEnergy =
            beamEnergy;

        // ----------------------------------------------------
        // Central 2x2 energy fraction
        //
        // Block layout:
        //
        // 0  1  2  3
        // 4  5  6  7
        // 8  9 10 11
        // 12 13 14 15
        //
        // Central blocks = 5,6,9,10
        // ----------------------------------------------------

        const G4double coreEnergy =
            blockEnergy[5] +
            blockEnergy[6] +
            blockEnergy[9] +
            blockEnergy[10];

        G4double coreFraction = 0.0;

        if (totalEnergy > 0.0)
        {
            coreFraction =
                coreEnergy / totalEnergy;
        }

        // ----------------------------------------------------
        // Lateral shower width
        // ----------------------------------------------------

        G4double weightedR2 = 0.0;

        for (G4int row = 0; row < 4; ++row)
        {
            for (G4int col = 0; col < 4; ++col)
            {
                const G4int index =
                    row * 4 + col;

                const G4double x =
                    (static_cast<G4double>(col) - 1.5)
                    * 10.0 * cm;

                const G4double y =
                    (static_cast<G4double>(row) - 1.5)
                    * 10.0 * cm;

                const G4double r2 =
                    x * x + y * y;

                weightedR2 +=
                    blockEnergy[index] * r2;
            }
        }

        G4double showerWidth = 0.0;

        if (totalEnergy > 0.0)
        {
            showerWidth =
                std::sqrt(
                    weightedR2 / totalEnergy
                );
        }

        // ----------------------------------------------------
        // Write EVERY event.
        //
        // No filtering happens here.
        // ----------------------------------------------------

        
        CsvOutput::Instance().WriteEvent(
            G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID(),
            event->GetEventID(),
            nominalBeamEnergy,
            material,
            absorberThickness,
            absorberX0,
            beamEnergy,
            beamMomentum,
            vertexX,
            vertexY,
            blockEnergy,
            totalEnergy,
            coreFraction,
            showerWidth
        );
    }

private:

    const DetectorConstruction* fDetector = nullptr;
};


// ============================================================
// ActionInitialisation
// ============================================================

void ActionInitialisation::BuildForMaster() const
{
    SetUserAction(
        new RunAction()
    );
}


void ActionInitialisation::Build() const
{
    auto* detector =
        static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()
                ->GetUserDetectorConstruction()
        );

    SetUserAction(
        new PrimaryGeneratorAction()
    );

    SetUserAction(
        new RunAction()
    );

    SetUserAction(
        new EventAction(detector)
    );
}
