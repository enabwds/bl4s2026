// ============================================================
//  ActionInitialisation.cc
//  Defines PrimaryGeneratorAction, RunAction, EventAction,
//  and NoiseFilter — all in one file for student readability.
// ============================================================

#include "ActionInitialisation.hh"
#include "CalorimeterSD.hh"
#include "DetectorConstruction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mutex>

// ============================================================
//  CSV OUTPUT
// ============================================================
namespace CsvOutput {
    std::ofstream file;
    std::mutex    mtx;

    void open(const std::string& fname) {
        if (file.is_open()) return;   // already open from a previous run — don't truncate
        bool needHeader = false;
        {
            std::ifstream probe(fname, std::ios::ate);
            needHeader = (!probe.is_open() || probe.tellg() == 0);
        }
        file.open(fname, std::ios::app);   // APPEND — never truncates existing data
        if (needHeader) {
            file << "event_id,material,absorber_thickness_mm,absorber_x0,"
                 << "beam_energy_GeV,beam_px_MeVc,beam_py_MeVc,beam_pz_MeVc,"
                 << "vertex_x_mm,vertex_y_mm,";
            for (int i = 0; i < 16; ++i)
                file << "E_block_" << i << "_GeV,";
            file << "Etotal_GeV,CoreFraction,ShowerWidth_cm\n";
        }
    }

    void write(long eventID,
               const std::string& material,
               double thickMM, double thickX0,
               double beamGeV,
               double px, double py, double pz,
               double vx, double vy,
               const std::array<double,16>& E,
               double Etot, double core, double width)
    {
        std::lock_guard<std::mutex> lock(mtx);
        file << eventID << ","
             << material << ","
             << std::fixed << std::setprecision(4)
             << thickMM  << ","
             << thickX0  << ","
             << beamGeV  << ","
             << px << "," << py << "," << pz << ","
             << vx << "," << vy << ",";
        for (double e : E) file << e << ",";
        file << Etot << "," << core << "," << width << "\n";
    }

    void close() { file.close(); }
}

// ============================================================
//  PRIMARY GENERATOR
//
//  Models the real BL4S PS beamline (T9/H4):
//
//  Beam species  : electrons (e-), >90% purity below 3 GeV
//                  (Cherenkov detectors confirm purity in simulation
//                   — here we just shoot pure electrons)
//
//  Energies      : 1, 2, 4 GeV (set via /gun/energy macro command)
//
//  Momentum spread: ±0.15 GeV/c, Gaussian (fixed by beamline optics)
//                   σ_p = 0.15/2.35 GeV/c (FWHM → sigma conversion)
//
//  Beam spot     : 2 cm diameter circular cross-section at focal point
//                  Modelled as uniform disk (conservative)
//
//  Divergence    : ~1 mrad (unconfirmed — check with BL4S team)
//                  Modelled as Gaussian angular spread σ = 0.5 mrad
//
//  Particles/spill: 10^4 to 10^5 (tunable via collimator)
//  DAQ rate limit : 3000 particles/second
//  Spill duration : ~400 ms  →  max ~1200 particles per spill at DAQ limit
// ============================================================
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction()
    {
        fGun = new G4ParticleGun(1);
        auto* table = G4ParticleTable::GetParticleTable();
        fGun->SetParticleDefinition(table->FindParticle("e-"));
        fGun->SetParticleEnergy(2.0 * GeV);
        fGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
        fGun->SetParticlePosition(G4ThreeVector(0, 0, -190.*cm));
    }

    ~PrimaryGeneratorAction() override { delete fGun; }

    void GeneratePrimaries(G4Event* event) override
    {
        // ── Nominal beam energy ───────────────────────────────────────
        G4double E0 = fGun->GetParticleEnergy();
        G4double m  = fGun->GetParticleDefinition()->GetPDGMass();
        G4double p0 = std::sqrt(E0*(E0 + 2*m));  // nominal momentum

        // ── Momentum spread: σ_p = 0.15 GeV/c / 2.355 (FWHM→σ) ─────
        // The beamline optics fix Δp/p ~ ±0.15 GeV/c FWHM
        G4double sigmaP = (0.15*GeV) / 2.355;
        G4double pz     = G4RandGauss::shoot(p0, sigmaP);
        pz = std::max(pz, 0.1*GeV);  // prevent negative momentum

        // ── Beam divergence: σ_angle = 0.5 mrad ─────────────────────
        // ~1 mrad total divergence (unconfirmed — update when confirmed)
        G4double sigmaTh = 0.5e-3;  // rad
        G4double theta_x = G4RandGauss::shoot(0., sigmaTh);
        G4double theta_y = G4RandGauss::shoot(0., sigmaTh);

        G4double ptot = pz / std::cos(std::sqrt(theta_x*theta_x + theta_y*theta_y));
        G4ThreeVector dir(std::sin(theta_x), std::sin(theta_y),
                          std::cos(std::sqrt(theta_x*theta_x + theta_y*theta_y)));
        dir = dir.unit();

        // ── Beam spot: uniform disk, 2 cm diameter ────────────────────
        G4double r   = 1.0*cm * std::sqrt(G4UniformRand());  // uniform in disk
        G4double phi = 2.*pi * G4UniformRand();
        G4double vx  = r * std::cos(phi);
        G4double vy  = r * std::sin(phi);

        // ── Kinetic energy from smeared momentum ──────────────────────
        G4double Esmeared = std::sqrt(ptot*ptot + m*m) - m;

        fGun->SetParticleEnergy(Esmeared);
        fGun->SetParticleMomentumDirection(dir);
        fGun->SetParticlePosition(G4ThreeVector(vx, vy, -190.*cm));
        fGun->GeneratePrimaryVertex(event);

        // Store vertex and momentum for CSV output
        fLastVx = vx/mm;
        fLastVy = vy/mm;
        fLastPx = ptot * dir.x() / (MeV);
        fLastPy = ptot * dir.y() / (MeV);
        fLastPz = ptot * dir.z() / (MeV);
    }

    // Accessors for EventAction
    double GetLastVx() const { return fLastVx; }
    double GetLastVy() const { return fLastVy; }
    double GetLastPx() const { return fLastPx; }
    double GetLastPy() const { return fLastPy; }
    double GetLastPz() const { return fLastPz; }

private:
    G4ParticleGun* fGun;
    double fLastVx=0, fLastVy=0;
    double fLastPx=0, fLastPy=0, fLastPz=0;
};

// ============================================================
//  RUN ACTION
// ============================================================
class RunAction : public G4UserRunAction
{
public:
    RunAction() = default;

    void BeginOfRunAction(const G4Run*) override {
        if (isMaster) {
            CsvOutput::open("shower_output.csv");
            G4cout << "Output: shower_output.csv\n";
        }
    }

    void EndOfRunAction(const G4Run* run) override {
        if (isMaster) {
            CsvOutput::file.flush();   // flush to disk; keep file open for next run
            G4cout << "Run complete. Events: " << run->GetNumberOfEvent() << "\n";
        }
    }
};

// ============================================================
//  NOISE FILTER
//
//  Three-level noise rejection matching realistic beamline conditions:
//
//  Cut 1 — Sub-threshold block deposits (in CalorimeterSD.cc):
//    Steps below 0.5 MeV dropped before accumulation.
//    Mimics PMT discriminator threshold on lead-glass blocks.
//
//  Cut 2 — Beam miss (event-level):
//    Etotal < 5% of beam energy → beam missed calorimeter.
//    In real experiment: DWC tracks would show no hit in calo acceptance.
//
//  Cut 3 — Statistical outlier (event-level, Welford algorithm):
//    Events > 3.5σ from running mean are anomalous secondaries
//    (e.g. back-scatter, δ-rays reaching detector from upstream).
//    Welford online algorithm: O(1) memory, numerically stable.
// ============================================================
struct NoiseFilter
{
    static constexpr double kMinFraction  = 0.05;
    static constexpr double kOutlierSigma = 3.5;

    long   n    = 0;
    double mean = 0.;
    double M2   = 0.;

    void Reset() {
        n    = 0;
        mean = 0.;
        M2   = 0.;
    }

    void update(double x) {
        ++n;
        double d1 = x - mean;
        mean += d1 / n;
        M2   += d1 * (x - mean);
    }

    double stddev() const { return (n > 1) ? std::sqrt(M2/(n-1)) : 0.; }

    bool accept(double Etot, double beamGeV) {
        if (Etot < kMinFraction * beamGeV) return false;
        if (n >= 10) {
            double s = stddev();
            if (s > 0. && std::abs(Etot - mean) > kOutlierSigma * s)
                return false;
        }
        update(Etot);
        return true;
    }
};

// ============================================================
//  EVENT ACTION
// ============================================================
class EventAction : public G4UserEventAction
{
public:
    EventAction() = default;
    ~EventAction() override = default;

    void EndOfEventAction(const G4Event* event) override
    {
        G4Event* ev = const_cast<G4Event*>(event);

        // ── Hit collection ───────────────────────────────────────────
        auto* hce = ev->GetHCofThisEvent();
        if (!hce) return;
        int hcID = G4SDManager::GetSDMpointer()->GetCollectionID("HitsCollection");
        auto* hc = dynamic_cast<CaloHitsCollection*>(hce->GetHC(hcID));
        if (!hc) return;

        // ── Block energies ───────────────────────────────────────────
        std::array<double,16> E = {};
        for (int i = 0; i < 16; ++i)
            E[i] = (*hc)[i]->GetEdep() / GeV;

        double Etot = 0.;
        for (double e : E) Etot += e;

        // ── Beam metadata ────────────────────────────────────────────
        double beamGeV = ev->GetPrimaryVertex()
                           ->GetPrimary()->GetKineticEnergy() / GeV;

        // ── Detector metadata ────────────────────────────────────────
        // Must be fetched BEFORE the noise filter so we can detect
        // configuration changes and reset the filter's running statistics.
        auto* det = static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        std::string material    = det->GetAbsorberMaterial()->GetName();
        double      thickMM     = det->GetAbsorberThickness();
        double      thickX0     = (det->GetAbsorberThickness()*mm) /
                                   det->GetAbsorberMaterial()->GetRadlen();

        // ── Reset noise filter on configuration change ────────────────
        // When material or thickness changes, the expected Etot distribution
        // shifts. The Welford running mean from the previous configuration
        // would incorrectly flag good events as outliers for ~50-100 events
        // at the start of each new configuration. Reset to avoid this.
        if (material != fLastMaterial ||
            std::abs(thickMM - fLastThickness) > 0.01) {
            fFilter.Reset();
            fLastMaterial  = material;
            fLastThickness = thickMM;
        }

        // ── Noise filter ─────────────────────────────────────────────
        if (!fFilter.accept(Etot, beamGeV)) { ++fRejected; return; }

        // ── Observables ──────────────────────────────────────────────

        // Core fraction: central 2×2 (blocks 5,6,9,10 in 4×4 grid)
        double Ecore        = E[5] + E[6] + E[9] + E[10];
        double coreFraction = (Etot > 0.) ? Ecore / Etot : 0.;

        // Lateral width: 2D energy-weighted radial RMS from array centre
        // Block centres (cm): col/row 0=-15, 1=-5, 2=+5, 3=+15
        // Using radial RMS from origin (standard Moliere radius proxy used
        // in real calorimeter analyses). The beam is well-centred so the
        // mean position is near zero — no mean subtraction needed.
        // This captures the full 2D lateral profile, not just the X axis.
        static const double xCentre[4] = {-15., -5.,  5., 15.};
        static const double yCentre[4] = {-15., -5.,  5., 15.};
        double sumW=0., sumWr2=0.;
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                double w = E[row*4 + col];
                double x = xCentre[col];
                double y = yCentre[row];
                sumW   += w;
                sumWr2 += w*(x*x + y*y);
            }
        }
        double width = (sumW > 0.) ? std::sqrt(sumWr2 / sumW) : 0.;

        // ── Beam vertex/momentum (from generator) ────────────────────
        // Access via primary vertex directly
        auto* vtx = ev->GetPrimaryVertex();
        double vx = vtx->GetX0()/mm;
        double vy = vtx->GetY0()/mm;
        auto*  p  = vtx->GetPrimary();
        double px = p->GetPx()/(MeV);
        double py = p->GetPy()/(MeV);
        double pz = p->GetPz()/(MeV);

        // ── Write CSV ────────────────────────────────────────────────
        CsvOutput::write(ev->GetEventID(),
                         material, thickMM, thickX0, beamGeV,
                         px, py, pz, vx, vy,
                         E, Etot, coreFraction, width);
    }

    long GetRejected() const { return fRejected; }

private:
    NoiseFilter fFilter;
    long        fRejected      = 0;
    std::string fLastMaterial  = "";
    double      fLastThickness = -1.;
};

// ── ActionInitialisation ──────────────────────────────────────────────────────
void ActionInitialisation::BuildForMaster() const { SetUserAction(new RunAction()); }

void ActionInitialisation::Build() const {
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new RunAction());
    SetUserAction(new EventAction());
}
