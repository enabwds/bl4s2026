#include "ActionInitialisation.hh"
#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"

#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4GenericMessenger.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4ios.hh"

#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>


// ============================================================
// Simulation configuration
// ============================================================

namespace SimulationConfig
{
    constexpr const char* kDatasetSchemaVersion =
        "BremScope_EventSchema_v3";

    constexpr const char* kSimulationConfiguration =
        "4x4_lead_glass_calorimeter";

    // --------------------------------------------------------
    // Beam configuration
    // --------------------------------------------------------

    constexpr G4double kBeamMomentumSpreadFWHM =
        0.15 * GeV;

    constexpr G4double kBeamAngularSigma =
        0.5 * mrad;

    constexpr G4double kBeamSpotRadius =
        1.0 * cm;

    constexpr G4double kBeamVertexZ =
        -191.0 * cm;

    // --------------------------------------------------------
    // Calorimeter configuration
    // --------------------------------------------------------

    constexpr G4int kNumberOfCalorimeterBlocks =
        16;

    constexpr G4int kCalorimeterRows =
        4;

    constexpr G4int kCalorimeterColumns =
        4;

    constexpr G4double kBlockPitch =
        10.0 * cm;

    constexpr G4double kBlockSize =
        10.0 * cm;

    constexpr G4double kCalorimeterDepth =
        37.0 * cm;

    // --------------------------------------------------------
    // Radial bins
    // --------------------------------------------------------

    constexpr G4double kRadialBin1 =
        5.0 * cm;

    constexpr G4double kRadialBin2 =
        10.0 * cm;

    constexpr G4double kRadialBin3 =
        15.0 * cm;

    constexpr G4double kRadialBin4 =
        20.0 * cm;
}


// ============================================================
// Per-event GEANT4 truth / accounting data
// ============================================================

struct TruthEventData
{
    G4double absorberEdep = 0.0;
    G4double otherDetectorEdep = 0.0;

    G4double calorimeterEscapeEnergy = 0.0;
    G4double worldEscapeEnergy = 0.0;

    std::unordered_set<G4int> electronTracks;
    std::unordered_set<G4int> positronTracks;
    std::unordered_set<G4int> gammaTracks;
    std::unordered_set<G4int> otherTracks;

    void Reset()
    {
        absorberEdep = 0.0;
        otherDetectorEdep = 0.0;

        calorimeterEscapeEnergy = 0.0;
        worldEscapeEnergy = 0.0;

        electronTracks.clear();
        positronTracks.clear();
        gammaTracks.clear();
        otherTracks.clear();
    }

    G4int GetElectronCount() const
    {
        return static_cast<G4int>(
            electronTracks.size()
        );
    }

    G4int GetPositronCount() const
    {
        return static_cast<G4int>(
            positronTracks.size()
        );
    }

    G4int GetGammaCount() const
    {
        return static_cast<G4int>(
            gammaTracks.size()
        );
    }

    G4int GetOtherCount() const
    {
        return static_cast<G4int>(
            otherTracks.size()
        );
    }

    G4int GetTotalTrackCount() const
    {
        return
            GetElectronCount()
            + GetPositronCount()
            + GetGammaCount()
            + GetOtherCount();
    }
};


// ============================================================
// Primary generator
// ============================================================

class PrimaryGeneratorAction
    : public G4VUserPrimaryGeneratorAction
{
public:

    PrimaryGeneratorAction()
    {
        fParticleGun =
            new G4ParticleGun(1);

        auto* electron =
            G4ParticleTable::GetParticleTable()
                ->FindParticle("e-");

        fParticleGun->SetParticleDefinition(
            electron
        );

        fNominalBeamEnergy =
            2.0 * GeV;

        fMessenger =
            new G4GenericMessenger(
                this,
                "/gun/",
                "Primary beam controls"
            );

        fMessenger->DeclarePropertyWithUnit(
            "setNominalEnergy",
            "GeV",
            fNominalBeamEnergy,
            "Set nominal beam energy"
        );

        fParticleGun->SetParticlePosition(
            G4ThreeVector(
                0.0,
                0.0,
                SimulationConfig::kBeamVertexZ
            )
        );
    }

    ~PrimaryGeneratorAction() override
    {
        delete fParticleGun;
        delete fMessenger;
    }

    void GeneratePrimaries(G4Event* event) override
    {
        const G4double nominalEnergy =
            fNominalBeamEnergy;

        const G4double sigmaP =
            SimulationConfig::kBeamMomentumSpreadFWHM
            / 2.354820045;

        G4double momentum =
            G4RandGauss::shoot(
                nominalEnergy,
                sigmaP
            );

        momentum =
            std::max(
                momentum,
                0.1 * GeV
            );

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

        direction =
            direction.unit();

        const G4double radius =
            SimulationConfig::kBeamSpotRadius
            * std::sqrt(G4UniformRand());

        const G4double phi =
            CLHEP::twopi
            * G4UniformRand();

        const G4double x =
            radius * std::cos(phi);

        const G4double y =
            radius * std::sin(phi);

        fParticleGun->SetParticleMomentum(
            momentum
        );

        fParticleGun->SetParticleMomentumDirection(
            direction
        );

        fParticleGun->SetParticlePosition(
            G4ThreeVector(
                x,
                y,
                SimulationConfig::kBeamVertexZ
            )
        );

        fParticleGun->GeneratePrimaryVertex(
            event
        );
    }

    G4double GetNominalBeamEnergy() const
    {
        return fNominalBeamEnergy;
    }

private:

    G4ParticleGun* fParticleGun =
        nullptr;

    G4GenericMessenger* fMessenger =
        nullptr;

    G4double fNominalBeamEnergy =
        2.0 * GeV;
};


// ============================================================
// CSV output
// ============================================================

class CsvOutput
{
public:

    static CsvOutput& Instance()
    {
        static CsvOutput instance;
        return instance;
    }

    void Open(
        const std::string& filename
    )
    {
        std::lock_guard<std::mutex> lock(
            fMutex
        );

        if (fFile.is_open())
            return;

        fFile.open(
            filename,
            std::ios::out | std::ios::trunc
        );

        if (!fFile.is_open())
        {
            G4Exception(
                "CsvOutput::Open",
                "Output001",
                FatalException,
                "Could not open shower output file."
            );

            return;
        }

        fFile
            << std::setprecision(10);

        WriteHeader();

        G4cout
            << "Opened CSV: "
            << filename
            << G4endl;
    }

    void Close()
    {
        std::lock_guard<std::mutex> lock(
            fMutex
        );

        if (fFile.is_open())
        {
            fFile.flush();
            fFile.close();
        }
    }

    void WriteEvent(
        G4int runID,
        G4int eventID,
        const std::string& datasetSchemaVersion,
        const std::string& simulationConfiguration,
        G4double nominalBeamEnergy,
        G4double beamEnergy,
        const G4ThreeVector& beamMomentum,
        G4double beamThetaX,
        G4double beamThetaY,
        G4double vertexX,
        G4double vertexY,
        G4double vertexZ,
        const std::string& material,
        G4double absorberThickness,
        G4double absorberX0,
        const std::vector<G4double>& blockEnergy,
        G4double totalEnergy,
        G4double coreFraction,
        G4double showerCentroidX,
        G4double showerCentroidY,
        G4double showerSigmaX,
        G4double showerSigmaY,
        G4double showerWidth,
        G4double radialEnergy0to5,
        G4double radialEnergy5to10,
        G4double radialEnergy10to15,
        G4double radialEnergy15to20,
        G4double radialEnergyAbove20,
        G4double maxBlockEnergy,
        G4double maxBlockFraction,
        G4int maxBlockID,
        const TruthEventData& truth
    )
    {
        std::lock_guard<std::mutex> lock(
            fMutex
        );

        if (!fFile.is_open())
            return;

        fFile
            << runID << ','
            << eventID << ','
            << datasetSchemaVersion << ','
            << simulationConfiguration << ',';

        fFile
            << nominalBeamEnergy / GeV << ','
            << beamEnergy / GeV << ','
            << beamMomentum.x() / MeV << ','
            << beamMomentum.y() / MeV << ','
            << beamMomentum.z() / MeV << ','
            << beamThetaX / mrad << ','
            << beamThetaY / mrad << ',';

        fFile
            << vertexX / mm << ','
            << vertexY / mm << ','
            << vertexZ / mm << ',';

        fFile
            << SimulationConfig::kBeamMomentumSpreadFWHM / GeV
            << ','
            << SimulationConfig::kBeamAngularSigma / mrad
            << ','
            << SimulationConfig::kBeamSpotRadius / cm
            << ',';

        fFile
            << SimulationConfig::kNumberOfCalorimeterBlocks
            << ','
            << SimulationConfig::kCalorimeterRows
            << ','
            << SimulationConfig::kCalorimeterColumns
            << ','
            << SimulationConfig::kBlockPitch / cm
            << ','
            << SimulationConfig::kBlockSize / cm
            << ','
            << SimulationConfig::kCalorimeterDepth / cm
            << ',';

        fFile
            << material << ','
            << absorberThickness / mm << ','
            << absorberX0 / mm << ',';

        for (G4double energy : blockEnergy)
        {
            fFile
                << energy / GeV
                << ',';
        }

        for (G4double energy : blockEnergy)
        {
            const G4double fraction =
                totalEnergy > 0.0
                    ? energy / totalEnergy
                    : 0.0;

            fFile
                << fraction
                << ',';
        }

        fFile
            << totalEnergy / GeV << ',';

        fFile
            << coreFraction << ','
            << showerCentroidX / cm << ','
            << showerCentroidY / cm << ','
            << showerSigmaX / cm << ','
            << showerSigmaY / cm << ','
            << showerWidth / cm << ',';

        fFile
            << radialEnergy0to5 / GeV << ','
            << radialEnergy5to10 / GeV << ','
            << radialEnergy10to15 / GeV << ','
            << radialEnergy15to20 / GeV << ','
            << radialEnergyAbove20 / GeV << ',';

        const G4double radialFraction0to5 =
            totalEnergy > 0.0
                ? radialEnergy0to5 / totalEnergy
                : 0.0;

        const G4double radialFraction5to10 =
            totalEnergy > 0.0
                ? radialEnergy5to10 / totalEnergy
                : 0.0;

        const G4double radialFraction10to15 =
            totalEnergy > 0.0
                ? radialEnergy10to15 / totalEnergy
                : 0.0;

        const G4double radialFraction15to20 =
            totalEnergy > 0.0
                ? radialEnergy15to20 / totalEnergy
                : 0.0;

        const G4double radialFractionAbove20 =
            totalEnergy > 0.0
                ? radialEnergyAbove20 / totalEnergy
                : 0.0;

        fFile
            << radialFraction0to5 << ','
            << radialFraction5to10 << ','
            << radialFraction10to15 << ','
            << radialFraction15to20 << ','
            << radialFractionAbove20 << ',';

        fFile
            << maxBlockEnergy / GeV << ','
            << maxBlockFraction << ','
            << maxBlockID << ',';

        fFile
            << truth.absorberEdep / GeV << ','
            << truth.otherDetectorEdep / GeV << ','
            << truth.calorimeterEscapeEnergy / GeV << ','
            << truth.worldEscapeEnergy / GeV << ',';

        fFile
            << truth.GetElectronCount() << ','
            << truth.GetPositronCount() << ','
            << truth.GetGammaCount() << ','
            << truth.GetOtherCount() << ','
            << truth.GetTotalTrackCount() << ',';

        const G4double accountedEnergy =
            totalEnergy
            + truth.absorberEdep
            + truth.otherDetectorEdep
            + truth.worldEscapeEnergy;

        const G4double residualEnergy =
            beamEnergy - accountedEnergy;

        fFile
            << accountedEnergy / GeV << ','
            << residualEnergy / GeV
            << '\n';
    }

private:

    CsvOutput() = default;

    ~CsvOutput()
    {
        Close();
    }

    void WriteHeader()
    {
        fFile
            << "run_id,"
            << "event_id,"
            << "dataset_schema_version,"
            << "simulation_configuration,";

        fFile
            << "nominal_beam_energy_GeV,"
            << "beam_energy_GeV,"
            << "beam_px_MeVc,"
            << "beam_py_MeVc,"
            << "beam_pz_MeVc,"
            << "beam_theta_x_mrad,"
            << "beam_theta_y_mrad,";

        fFile
            << "vertex_x_mm,"
            << "vertex_y_mm,"
            << "vertex_z_mm,";

        fFile
            << "beam_momentum_spread_FWHM_GeV,"
            << "beam_angular_sigma_mrad,"
            << "beam_spot_radius_cm,";

        fFile
            << "calorimeter_nblocks,"
            << "calorimeter_nrows,"
            << "calorimeter_ncolumns,"
            << "block_pitch_cm,"
            << "block_size_cm,"
            << "calorimeter_depth_cm,";

        fFile
            << "material,"
            << "absorber_thickness_mm,"
            << "absorber_x0_mm";

        for (
            G4int i = 0;
            i < SimulationConfig::kNumberOfCalorimeterBlocks;
            ++i
        )
        {
            fFile
                << ",E_block_"
                << i
                << "_GeV";
        }

        for (
            G4int i = 0;
            i < SimulationConfig::kNumberOfCalorimeterBlocks;
            ++i
        )
        {
            fFile
                << ",F_block_"
                << i;
        }

        fFile
            << ",Etotal_GeV";

        fFile
            << ",CoreFraction"
            << ",ShowerCentroidX_cm"
            << ",ShowerCentroidY_cm"
            << ",ShowerSigmaX_cm"
            << ",ShowerSigmaY_cm"
            << ",ShowerWidth_cm";

        fFile
            << ",E_radial_0to5cm_GeV"
            << ",E_radial_5to10cm_GeV"
            << ",E_radial_10to15cm_GeV"
            << ",E_radial_15to20cm_GeV"
            << ",E_radial_above20cm_GeV";

        fFile
            << ",F_radial_0to5cm"
            << ",F_radial_5to10cm"
            << ",F_radial_10to15cm"
            << ",F_radial_15to20cm"
            << ",F_radial_above20cm";

        fFile
            << ",MaxBlockEnergy_GeV"
            << ",MaxBlockFraction"
            << ",MaxBlockID";

        fFile
            << ",E_absorber_GeV"
            << ",E_other_detector_GeV"
            << ",E_escape_calorimeter_GeV"
            << ",E_escape_world_GeV";

        fFile
            << ",N_electron_tracks"
            << ",N_positron_tracks"
            << ",N_gamma_tracks"
            << ",N_other_tracks"
            << ",N_total_tracks";

        fFile
            << ",E_kinetic_accounted_GeV"
            << ",E_kinetic_residual_GeV"
            << '\n';

        fFile.flush();
    }

    std::ofstream fFile;
    std::mutex fMutex;
};


// ============================================================
// Run action
// ============================================================

class RunAction
    : public G4UserRunAction
{
public:

    ~RunAction() override
    {
        if (
            !G4Threading::IsMultithreadedApplication()
            ||
            G4Threading::IsMasterThread()
        )
        {
            CsvOutput::Instance().Close();
        }
    }

    void BeginOfRunAction(
        const G4Run* run
    ) override
    {
        if (!run)
            return;

        const bool isMT =
            G4Threading::IsMultithreadedApplication();

        const bool isMaster =
            G4Threading::IsMasterThread();

        if (isMT && isMaster)
            return;

        const std::string filename =
            "shower_output_run_"
            + std::to_string(
                run->GetRunID()
            )
            + ".csv";

        CsvOutput::Instance()
            .Open(filename);
    }

    void EndOfRunAction(
        const G4Run*
    ) override
    {
        const bool isMT =
            G4Threading::IsMultithreadedApplication();

        const bool isMaster =
            G4Threading::IsMasterThread();

        if (!isMT || isMaster)
        {
            CsvOutput::Instance().Close();
        }
    }
};


// ============================================================
// Tracking action
// ============================================================

class TrackingAction
    : public G4UserTrackingAction
{
public:

    explicit TrackingAction(
        const std::shared_ptr<TruthEventData>& truth
    )
        : fTruth(truth)
    {
    }

    void PreUserTrackingAction(
        const G4Track* track
    ) override
    {
        if (!track || !fTruth)
            return;

        const G4int trackID =
            track->GetTrackID();

        if (trackID < 0)
            return;

        const auto* particle =
            track->GetParticleDefinition();

        if (!particle)
            return;

        const G4String particleName =
            particle->GetParticleName();

        if (particleName == "e-")
        {
            fTruth->electronTracks.insert(
                trackID
            );
        }
        else if (particleName == "e+")
        {
            fTruth->positronTracks.insert(
                trackID
            );
        }
        else if (particleName == "gamma")
        {
            fTruth->gammaTracks.insert(
                trackID
            );
        }
        else
        {
            fTruth->otherTracks.insert(
                trackID
            );
        }
    }

private:

    std::shared_ptr<TruthEventData> fTruth;
};


// ============================================================
// Stepping action
// ============================================================

class SteppingAction
    : public G4UserSteppingAction
{
public:

    explicit SteppingAction(
        const std::shared_ptr<TruthEventData>& truth
    )
        : fTruth(truth)
    {
    }

    void UserSteppingAction(
        const G4Step* step
    ) override
    {
        if (!step || !fTruth)
            return;

        const auto* prePoint =
            step->GetPreStepPoint();

        const auto* postPoint =
            step->GetPostStepPoint();

        const auto* track =
            step->GetTrack();

        if (!prePoint || !postPoint || !track)
            return;

        const G4double edep =
            step->GetTotalEnergyDeposit();

        const auto* preVolume =
            prePoint->GetPhysicalVolume();

        const auto* postVolume =
            postPoint->GetPhysicalVolume();

        const G4String preVolumeName =
            preVolume
                ? preVolume->GetLogicalVolume()->GetName()
                : "";

        const G4String postVolumeName =
            postVolume
                ? postVolume->GetLogicalVolume()->GetName()
                : "";

        const bool inAbsorber =
            preVolumeName == "Absorber";

        const bool inCalorimeter =
            preVolumeName.find("Block_") == 0;

        // ----------------------------------------------------
        // TEMPORARY ABSORBER DEBUG
        //
        // Print the first 20 steps occurring inside the
        // absorber. This is only for the sanity check and
        // should be removed afterwards.
        // ----------------------------------------------------
 
        // ----------------------------------------------------
        // Absorber energy deposition
        // ----------------------------------------------------

        if (edep > 0.0)
        {
            if (inAbsorber)
            {
                fTruth->absorberEdep +=
                    edep;
            }
            else if (!inCalorimeter)
            {
                fTruth->otherDetectorEdep +=
                    edep;
            }
        }

        // ----------------------------------------------------
        // Calorimeter escape
        // ----------------------------------------------------

        if (inCalorimeter)
        {
            const bool stillInCalorimeter =
                postVolume
                &&
                postVolumeName.find("Block_") == 0;

            if (!stillInCalorimeter)
            {
                const G4double kineticEnergy =
                    track->GetKineticEnergy();

                if (kineticEnergy > 0.0)
                {
                    fTruth->calorimeterEscapeEnergy +=
                        kineticEnergy;
                }
            }
        }

        // ----------------------------------------------------
        // World escape
        // ----------------------------------------------------

        if (
            postPoint->GetStepStatus()
            == fWorldBoundary
        )
        {
            const G4double kineticEnergy =
                track->GetKineticEnergy();

            if (kineticEnergy > 0.0)
            {
                fTruth->worldEscapeEnergy +=
                    kineticEnergy;
            }
        }
    }

private:

    std::shared_ptr<TruthEventData> fTruth;
};


// ============================================================
// Event action
// ============================================================

class EventAction
    : public G4UserEventAction
{
public:

    EventAction(
        const DetectorConstruction* detector,
        const PrimaryGeneratorAction* primaryGenerator,
        const std::shared_ptr<TruthEventData>& truth
    )
        :
        fDetector(detector),
        fPrimaryGenerator(primaryGenerator),
        fTruth(truth)
    {
    }

    void BeginOfEventAction(
        const G4Event*
    ) override
    {
        if (fTruth)
        {
            fTruth->Reset();
        }
    }

    void EndOfEventAction(
        const G4Event* event
    ) override
    {
        if (!event)
            return;

        if (!fTruth)
            return;

        auto* hce =
            event->GetHCofThisEvent();

        if (!hce)
            return;

        auto* sdManager =
            G4SDManager::GetSDMpointer();

        const G4int collectionID =
            sdManager->GetCollectionID(
                "CalorimeterSD/CaloHits"
            );

        if (collectionID < 0)
            return;

        auto* hits =
            static_cast<CaloHitsCollection*>(
                hce->GetHC(collectionID)
            );

        if (!hits)
            return;

        std::vector<G4double> blockEnergy(
            SimulationConfig::kNumberOfCalorimeterBlocks,
            0.0
        );

        G4double totalEnergy =
            0.0;

        const G4int numberToRead =
            std::min(
                static_cast<G4int>(
                    hits->entries()
                ),
                SimulationConfig::kNumberOfCalorimeterBlocks
            );

        for (
            G4int i = 0;
            i < numberToRead;
            ++i
        )
        {
            auto* hit =
                (*hits)[i];

            if (!hit)
                continue;

            const G4int blockID =
                hit->GetBlockID();

            if (
                blockID < 0
                ||
                blockID >=
                    SimulationConfig::kNumberOfCalorimeterBlocks
            )
            {
                continue;
            }

            blockEnergy[blockID] =
                hit->GetEdep();

            totalEnergy +=
                blockEnergy[blockID];
        }

        const auto* primaryVertex =
            event->GetPrimaryVertex();

        if (!primaryVertex)
            return;

        const auto* primary =
            primaryVertex->GetPrimary();

        if (!primary)
            return;

        const G4double beamEnergy =
            primary->GetKineticEnergy();

        const G4ThreeVector beamMomentum =
            primary->GetMomentum();

        const G4double vertexX =
            primaryVertex->GetX0();

        const G4double vertexY =
            primaryVertex->GetY0();

        const G4double vertexZ =
            primaryVertex->GetZ0();

        const G4double nominalBeamEnergy =
            fPrimaryGenerator
                ->GetNominalBeamEnergy();

        G4double beamThetaX =
            0.0;

        G4double beamThetaY =
            0.0;

        if (std::abs(beamMomentum.z()) > 0.0)
        {
            beamThetaX =
                std::atan2(
                    beamMomentum.x(),
                    beamMomentum.z()
                );

            beamThetaY =
                std::atan2(
                    beamMomentum.y(),
                    beamMomentum.z()
                );
        }

        const auto* absorberMaterial =
            fDetector->GetAbsorberMaterial();

        if (!absorberMaterial)
            return;

        const std::string material =
            absorberMaterial->GetName();

        const G4double absorberThickness =
            fDetector->GetAbsorberThickness();

        const G4double absorberX0 =
            absorberMaterial->GetRadlen();

        const G4double coreEnergy =
            blockEnergy[5]
            + blockEnergy[6]
            + blockEnergy[9]
            + blockEnergy[10];

        G4double coreFraction =
            0.0;

        if (totalEnergy > 0.0)
        {
            coreFraction =
                coreEnergy / totalEnergy;
        }

        G4double weightedX =
            0.0;

        G4double weightedY =
            0.0;

        for (
            G4int row = 0;
            row < SimulationConfig::kCalorimeterRows;
            ++row
        )
        {
            for (
                G4int col = 0;
                col < SimulationConfig::kCalorimeterColumns;
                ++col
            )
            {
                const G4int index =
                    row
                    * SimulationConfig::kCalorimeterColumns
                    + col;

                const G4double x =
                    (
                        static_cast<G4double>(col)
                        - 1.5
                    )
                    * SimulationConfig::kBlockPitch;

                const G4double y =
                    (
                        static_cast<G4double>(row)
                        - 1.5
                    )
                    * SimulationConfig::kBlockPitch;

                weightedX +=
                    blockEnergy[index] * x;

                weightedY +=
                    blockEnergy[index] * y;
            }
        }

        G4double showerCentroidX =
            0.0;

        G4double showerCentroidY =
            0.0;

        if (totalEnergy > 0.0)
        {
            showerCentroidX =
                weightedX / totalEnergy;

            showerCentroidY =
                weightedY / totalEnergy;
        }

        G4double weightedX2 =
            0.0;

        G4double weightedY2 =
            0.0;

        G4double weightedRadiusSquared =
            0.0;

        G4double radialEnergy0to5 =
            0.0;

        G4double radialEnergy5to10 =
            0.0;

        G4double radialEnergy10to15 =
            0.0;

        G4double radialEnergy15to20 =
            0.0;

        G4double radialEnergyAbove20 =
            0.0;

        G4double maxBlockEnergy =
            0.0;

        G4int maxBlockID =
            -1;

        for (
            G4int row = 0;
            row < SimulationConfig::kCalorimeterRows;
            ++row
        )
        {
            for (
                G4int col = 0;
                col < SimulationConfig::kCalorimeterColumns;
                ++col
            )
            {
                const G4int index =
                    row
                    * SimulationConfig::kCalorimeterColumns
                    + col;

                const G4double energy =
                    blockEnergy[index];

                const G4double x =
                    (
                        static_cast<G4double>(col)
                        - 1.5
                    )
                    * SimulationConfig::kBlockPitch;

                const G4double y =
                    (
                        static_cast<G4double>(row)
                        - 1.5
                    )
                    * SimulationConfig::kBlockPitch;

                const G4double dx =
                    x - showerCentroidX;

                const G4double dy =
                    y - showerCentroidY;

                const G4double radius =
                    std::sqrt(
                        dx * dx
                        + dy * dy
                    );

                weightedX2 +=
                    energy * dx * dx;

                weightedY2 +=
                    energy * dy * dy;

                weightedRadiusSquared +=
                    energy
                    * radius
                    * radius;

                if (
                    radius
                    < SimulationConfig::kRadialBin1
                )
                {
                    radialEnergy0to5 +=
                        energy;
                }
                else if (
                    radius
                    < SimulationConfig::kRadialBin2
                )
                {
                    radialEnergy5to10 +=
                        energy;
                }
                else if (
                    radius
                    < SimulationConfig::kRadialBin3
                )
                {
                    radialEnergy10to15 +=
                        energy;
                }
                else if (
                    radius
                    < SimulationConfig::kRadialBin4
                )
                {
                    radialEnergy15to20 +=
                        energy;
                }
                else
                {
                    radialEnergyAbove20 +=
                        energy;
                }

                if (energy > maxBlockEnergy)
                {
                    maxBlockEnergy =
                        energy;

                    maxBlockID =
                        index;
                }
            }
        }

        G4double showerSigmaX =
            0.0;

        G4double showerSigmaY =
            0.0;

        G4double showerWidth =
            0.0;

        if (totalEnergy > 0.0)
        {
            showerSigmaX =
                std::sqrt(
                    weightedX2
                    / totalEnergy
                );

            showerSigmaY =
                std::sqrt(
                    weightedY2
                    / totalEnergy
                );

            showerWidth =
                std::sqrt(
                    weightedRadiusSquared
                    / totalEnergy
                );
        }

        G4double maxBlockFraction =
            0.0;

        if (totalEnergy > 0.0)
        {
            maxBlockFraction =
                maxBlockEnergy
                / totalEnergy;
        }

        const G4Run* currentRun =
            G4RunManager::GetRunManager()
                ->GetCurrentRun();

        if (!currentRun)
            return;

        CsvOutput::Instance().WriteEvent(

            currentRun->GetRunID(),
            event->GetEventID(),
            SimulationConfig::kDatasetSchemaVersion,
            SimulationConfig::kSimulationConfiguration,

            nominalBeamEnergy,
            beamEnergy,
            beamMomentum,
            beamThetaX,
            beamThetaY,

            vertexX,
            vertexY,
            vertexZ,

            material,
            absorberThickness,
            absorberX0,

            blockEnergy,
            totalEnergy,

            coreFraction,
            showerCentroidX,
            showerCentroidY,
            showerSigmaX,
            showerSigmaY,
            showerWidth,

            radialEnergy0to5,
            radialEnergy5to10,
            radialEnergy10to15,
            radialEnergy15to20,
            radialEnergyAbove20,

            maxBlockEnergy,
            maxBlockFraction,
            maxBlockID,

            *fTruth
        );
    }

private:

    const DetectorConstruction* fDetector =
        nullptr;

    const PrimaryGeneratorAction*
        fPrimaryGenerator = nullptr;

    std::shared_ptr<TruthEventData> fTruth;
};


// ============================================================
// Action initialization
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

    auto truth =
        std::make_shared<TruthEventData>();

    auto* primaryGenerator =
        new PrimaryGeneratorAction();

    SetUserAction(
        primaryGenerator
    );

    SetUserAction(
        new RunAction()
    );

    SetUserAction(
        new EventAction(
            detector,
            primaryGenerator,
            truth
        )
    );

    SetUserAction(
        new SteppingAction(
            truth
        )
    );

    SetUserAction(
        new TrackingAction(
            truth
        )
    );
}
