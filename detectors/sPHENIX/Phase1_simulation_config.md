# Phase 1 Simulation Configuration

This table documents the current pp minbias truth-point/V0 study sample. Entries marked "to confirm" should be tightened before using this as final analysis documentation.

| Component | Information |
| --- | --- |
| Collision system | pp |
| Collision energy | sqrt(s) = 200 GeV |
| Generator | PYTHIA 8.312 through `PHPythia8`; HepMC written with `Fun4AllHepMCOutputManager` in `GeneratePythia8HepMC.C`. |
| Generator configuration | Detroit minimum-bias config: `$CALIBRATIONROOT/Generators/phpythia8_detroit_minBias.cfg`. Main settings: `Beams:idA = 2212`, `Beams:idB = 2212`, `Beams:eCM = 200`, `PDF:pSet = 17`, `MultipartonInteractions:bProfile = 2`, `MultipartonInteractions:ecmRef = 200`, `MultipartonInteractions:pT0Ref = 1.40`, `MultipartonInteractions:ecmPow = 0.135`, `MultipartonInteractions:coreRadius = 0.56`, `MultipartonInteractions:coreFraction = 0.78`, `ColourReconnection:range = 5.4`, `SoftQCD:inelastic = on`, `ParticleDecays:limitTau0 = on`, `ParticleDecays:tau0Max = 0.0000001`. PYTHIA initialization reports pp at 200 GeV with non-diffractive, single-diffractive, double-diffractive, and central-diffractive inelastic processes. |
| Event sample | HepMC generation: 100,000 events, 100 files x 1000 events/file (`pp_minbias_hepmc_100k`). G4 truth-point transport: 100,000 events, 10,000 DST files x 10 events/file (`pp_minbias_truthdst_100k_zero_vtx`). Generator-level truth V0 tree: 93,784 charged two-body decays in `truthV0Tree` (`K0S = 62,471`, `Lambda = 21,633`, `anti-Lambda = 9,680`). Reconstructed truth-point pair tree: 1,545,227 pair candidates in `pairTree`. |
| Detector transport | GEANT4 10.7.4 / geant4-10-07-patch-04. Physics list: `FTFP_BERT`. EvtGen decayer enabled; PHOTOS 3.64 appears in the run log. |
| Geometry | sPHENIX geometry from the local `Fun4All_G4_sPHENIX.C`/`G4Setup_sPHENIX.C` setup. Run log records `RUNNUMBER = 1`, `CDB_GLOBALTAG = MDC2`, `TIMESTAMP = 6`, world material `G4_AIR`, world shape `G4Tubs`. Stored run nodes include `PIPE`, `TPC`, `TPC_ENDCAP`, `MAGNET`, `BH_1`, forward/backward BH nodes, `RECO_TRACKING_GEOMETRY/TPCGEOMCONTAINER`, and `GEOMETRY_IO`. Field map CDB domain: `FIELDMAP_GAP`, file `sphenix3dbigmapxyz_gap_rebuild_v2.root`, scale factor 1. TPC FEE channel map domain: `TPC_FEE_CHANNEL_MAP`, file `TpcFeePadPlacementv4-ana494_2024p022-ana500.root`. Exact detector-geometry campaign tag beyond CDB `MDC2/timestamp 6`: to confirm. |
| Reconstruction | No standard TPC digitization, clustering, tracking, or vertexing in this truth-point sample. `PHG4TpcTruthPointBuilder` writes ideal TPC layer-intersection truth points to `G4HIT_TPC_TRUECLUSTER`; `PHG4TpcTruthPointTree` writes the flat truth-point tree. `TpcTruthEventTree` writes one event-level generator QA row per event to `truthEventTree`. `TpcTruthV0DecayTree` writes generator-level true V0 decays to `truthV0Tree`. `TpcV0CandidateTree` later fits grouped truth points with the current helix/PCA workflow and writes `pairTree`. |
| Event conditions | HepMC input vertices kept unchanged in G4 (`use_beam_vertex = false`). The generated HepMC files used here have zero primary vertex coordinates, so the effective primary-vertex spread is zero for this campaign. `Input::PILEUPRATE = 0`, so no pileup. Detector inefficiencies, detector hit smearing, standard TPC digitization, and pattern-recognition failures are not included in this ideal truth-point sample. |
| Randomness | HepMC/PYTHIA campaign used fixed per-job seeds: `seed = 100000 + process_id`; job 0 used seed 100000. `PHPythia8` derives/logs its internal seed from `PHRandomSeed` (job 0 log shows `PHPythia8 random seed: 434894410`). G4 campaign did not set a fixed `RANDOMSEED` flag in `Fun4All_G4_sPHENIX.C`; `PHRandomSeed::Verbosity(1)` logs seeds per job, apparently drawn from the default source. Example G4 job 0 first logged seed: `501660783`. For strict reproducibility, record all job logs or set deterministic G4 seeds in the campaign. |
| Software environment | The scripts used `sphenix_setup.sh new` plus local install `/sphenix/user/dcxchenxi/install`, so the resolved release depends on when each campaign ran. The HepMC generation campaign `pp_minbias_hepmc_100k` used `release_new/new.8` according to `log/pp_minbias_hepmc_100k/hepmc_0.out`. The DST-producing G4 truth-point campaign `pp_minbias_truthdst_100k_zero_vtx` used `release_new/new.10` according to `log/pp_minbias_truthdst_100k_zero_vtx/general_0.out`. As of June 21, 2026, `sphenix_setup.sh new` resolves `OFFLINE_MAIN` to `release_new/new.15`. The reproduction runners are now pinned to `sphenix_setup.sh -n new.15` rather than the moving `new` alias. |

To regenerate the Phase 1 sample with a consistent fixed release, start by regenerating the HepMC input:

```bash
cd /sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX

./submit_hepmc_generation_campaign.sh \
  pp_minbias_hepmc_100k_new15 \
  100000 \
  1000 \
  100000
```

This submits 100 Condor jobs, each generating 1000 pp minbias events. The expected outputs are:

```text
output/pp_minbias_hepmc_100k_new15/hepmc_files.list
output/pp_minbias_hepmc_100k_new15/completed/pp_minbias_*.hepmc
log/pp_minbias_hepmc_100k_new15/
```

After the jobs finish, check the file count and release log:

```bash
find output/pp_minbias_hepmc_100k_new15/completed \
  -maxdepth 1 -name 'pp_minbias_*.hepmc' | wc -l

rg -n "release_new/new.15|GeneratePythia8HepMC|PYTHIA" \
  log/pp_minbias_hepmc_100k_new15/hepmc_0.out
```

Useful campaign commands that produced this sample:

```bash
./submit_hepmc_generation_campaign.sh \
  pp_minbias_hepmc_100k \
  100000 \
  1000 \
  100000

./submit_general_campaign.sh \
  pp_minbias_truth_100k_zero_vtx \
  output/pp_minbias_hepmc_100k/hepmc_files.list \
  100000 \
  10 \
  1000 \
  false
```

Phase 1 generator-QA trees can be produced from the truth-point DST file list:

```bash
cd /sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX

root.exe -l -b -q 'Fun4All_Phase1Truth.C(100000, \
  "output/pp_minbias_truthdst_100k_zero_vtx/dst_files.list", \
  "output/phase1_truth_events_100k_zero_vtx.root", \
  "output/phase1_truth_v0_decays_100k_zero_vtx.root")'
```

For the full 100k sample, the Condor wrapper is preferred because the input is split into 10,000 small DST files:

```bash
cd /sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX

./submit_phase1_truth_campaign.sh \
  phase1_truth_100k_zero_vtx \
  output/pp_minbias_truthdst_100k_zero_vtx/dst_files.list \
  100 \
  10 \
  10000
```

This writes chunk outputs under:

```text
output/phase1_truth_100k_zero_vtx/completed/
```

Merge the event-level and V0-level outputs separately:

```bash
find output/phase1_truth_100k_zero_vtx/completed \
  -maxdepth 1 -name 'phase1_truth_events_*.root' \
  | sort > output/phase1_truth_100k_zero_vtx/event_files.list

find output/phase1_truth_100k_zero_vtx/completed \
  -maxdepth 1 -name 'phase1_truth_v0_decays_*.root' \
  | sort > output/phase1_truth_100k_zero_vtx/v0_files.list

hadd -f -j 8 \
  output/phase1_truth_events_100k_zero_vtx_merged.root \
  @output/phase1_truth_100k_zero_vtx/event_files.list

hadd -f -j 8 \
  output/phase1_truth_v0_decays_100k_zero_vtx_merged.root \
  @output/phase1_truth_100k_zero_vtx/v0_files.list
```

This writes:

| File | Tree | Purpose |
| --- | --- | --- |
| `phase1_truth_events_100k_zero_vtx.root` | `truthEventTree` | one row per event: primary vertex, charged multiplicity, event-level V0 yields |
| `phase1_truth_v0_decays_100k_zero_vtx.root` | `truthV0Tree` | one row per true charged two-body `K0S`, `Lambda`, or `anti-Lambda` decay |

The first Phase 1 reference plots can be produced with:

```bash
root -l -b -q 'drawPhase1GeneratorQA.C( \
  "output/phase1_truth_events_100k_zero_vtx_merged.root", \
  "output/phase1_truth_v0_decays_100k_zero_vtx_merged.root", \
  "output/phase1_generatorQA_100k_zero_vtx")'
```

The current event-level fiducial definitions stored in `truthEventTree` are:

| Quantity | Definition |
| --- | --- |
| `n_charged_primary_eta_pt` | generated/input primary charged particles with `abs(eta) < 1.1` and `pT > 0.2 GeV/c` |
| `n_charged_primary_no_daughters_eta_pt` | same, but excluding truth particles with daughters in the truth map |
| `n_*_charged_decay` | generated V0 parents with the expected charged two-body daughter pair |
| `n_*_primary_charged_decay` | same, but requiring the V0 parent track id to be positive in `PHG4TruthInfoContainer` |
| `n_*_fiducial` | charged-decay V0 parents with `abs(y) < 1.1` and `pT > 0 GeV/c` |

Important truth-record nuance: `PHG4TruthInfoContainer` is the post-G4 truth container, not a raw HepMC table. Some weak-decay V0 parents appear with negative track IDs even though they are physics decays, so `truthV0Tree` stores `parent_is_primary` and `parent_is_sphenix_primary` flags. Phase 1 plots should state whether they use all charged V0 decays, positive-track-ID primary V0s, or sPHENIX-primary V0s.
