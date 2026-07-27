#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include <GlobalVariables.C>

#include "G4Setup_sPHENIX.C"

#include <DisplayOn.C>
#include <G4_Mbd.C>
#include <G4_CaloTrigger.C>
#include <G4_Centrality.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_KFParticle.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>

#include <G4_User.C>
#include <QA.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4tpc/PHG4TpcTruthPointBuilder.h>
#include <g4tpc/PHG4TpcTruthPointTree.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <Rtypes.h>  // resolves R__LOAD_LIBRARY for clang-tidy
#include <TROOT.h>
#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4tpc.so)

// For HepMC Hijing
// try inputFile = /sphenix/sim/sim01/sphnxpro/sHijing_HepMC/sHijing_0-12fm.dat

int Fun4All_G4_sPHENIX(
    const int nEvents = 1,
    const std::string &inputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/pp_minbias_test.hepmc",
    const std::string &outputFile = "G4sPHENIX_TpcTruthPoints.root",
    const std::string &embed_input_file = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const int skip = 0,
    const std::string &outdir = ".",
    const bool useBeamVertex = false,
    const std::string &truthPointDstFile = "",
    const std::string &trkrHitSetDstFile = "",
    const bool useSimplePion = false,
    const int simplePionPdgId = 211,
    const double simplePionPt = 0.3,
    const double simplePionEta = 0.0,
    const double simplePionPhi = 0.0,
    const double simpleVertexZ = 0.0)
{
  if (useSimplePion &&
      ((simplePionPdgId != 211 && simplePionPdgId != -211) || simplePionPt <= 0.0))
  {
    std::cout << "Fun4All_G4_sPHENIX: invalid fixed-pion configuration" << std::endl;
    return 1;
  }

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  PHRandomSeed::Verbosity(1);
  CDBInterface::instance()->Verbosity(1);
  
  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  if (const char *seed_text = gSystem->Getenv("SIM_G4_RANDOM_SEED"))
  {
    try
    {
      const long seed = std::stol(seed_text);
      if (seed <= 0 || seed > std::numeric_limits<int>::max())
      {
        throw std::out_of_range("seed outside positive int range");
      }
      rc->set_IntFlag("RANDOMSEED", static_cast<int>(seed));
      std::cout << "Fun4All_G4_sPHENIX: fixed G4 RANDOMSEED=" << seed << std::endl;
    }
    catch (const std::exception &error)
    {
      std::cout << "Fun4All_G4_sPHENIX: invalid SIM_G4_RANDOM_SEED='"
                << seed_text << "': " << error.what() << std::endl;
      return 1;
    }
  }
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as seed
  // You can either set this to a random value using PHRandomSeed()
  // which will make all seeds identical (not sure what the point of
  // this would be:
  //  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // or set it to a fixed value so you can debug your code
  //  rc->set_IntFlag("RANDOMSEED", 12345);


  //===============
  // Input options
  //===============
  // verbosity setting (applies to all input managers)
  Input::VERBOSITY = 0;
  // First enable the input generators
  // Either:
  // read previously generated g4-hits files, in this case it opens a DST and skips
  // the simulations step completely. The G4Setup macro is only loaded to get information
  // about the number of layers used for the cell reco code
  //  Input::READHITS = true;
  INPUTREADHITS::filename[0] = inputFile;
  // if you use a filelist
  // INPUTREADHITS::listfile[0] = inputFile;
  // Or:
  // Use particle generator
  // And
  // Further choose to embed newly simulated events to a previous simulation. Not compatible with `readhits = true`
  // In case embedding into a production output, please double check your G4Setup_sPHENIX.C and G4_*.C consistent with those in the production macro folder
  // E.g. /sphenix/sim//sim01/production/2016-07-21/single_particle/spacal2d/
  //  Input::EMBED = true;
  INPUTEMBED::filename[0] = embed_input_file;
  // if you use a filelist
  //INPUTEMBED::listfile[0] = embed_input_file;

  Input::SIMPLE = useSimplePion;
  // Input::SIMPLE_NUMBER = 2; // if you need 2 of them
  // Input::SIMPLE_VERBOSITY = 1;

  // Enable this is emulating the nominal pp/pA/AA collision vertex distribution
  Input::BEAM_CONFIGURATION = Input::pp_COLLISION; // Input::AA_COLLISION (default), Input::pA_COLLISION, Input::pp_COLLISION

  //  Input::PYTHIA6 = true;

  // Input::PYTHIA8 = true;

  //  Input::GUN = true;
  //  Input::GUN_NUMBER = 3; // if you need 3 of them
  // Input::GUN_VERBOSITY = 1;

  // Input::COSMIC = true;

  //D0 generator
  //Input::DZERO = false;
  //Input::DZERO_VERBOSITY = 0;
  //Lambda_c generator //Not ready yet
  //Input::LAMBDAC = false;
  //Input::LAMBDAC_VERBOSITY = 0;
  // Upsilon generator
  //Input::UPSILON = true;
  //Input::UPSILON_NUMBER = 3; // if you need 3 of them
  //Input::UPSILON_VERBOSITY = 0;

  Input::HEPMC = !useSimplePion;
  if (Input::HEPMC)
  {
    INPUTHEPMC::filename = inputFile;
  }
  //-----------------
  // Hijing options (symmetrize hijing, add flow, add fermi motion)
  //-----------------
  //  INPUTHEPMC::HIJINGFLIP = true;
  //  INPUTHEPMC::FLOW = true;
  //  INPUTHEPMC::FLOW_VERBOSITY = 3;
  //  INPUTHEPMC::FERMIMOTION = true;


  // Event pile up simulation with collision rate in Hz MB collisions.
  //Input::PILEUPRATE = 50e3; // 50 kHz for AuAu
  //Input::PILEUPRATE = 3e6; // 3MHz for pp
  Input::PILEUPRATE = 0; // no pileup for this first TPC truth-point test

  // Enable this is emulating the nominal pp/pA/AA collision vertex distribution
  // for HepMC records (hijing, pythia8)
  //  Input::BEAM_CONFIGURATION = Input::AA_COLLISION; // for 2023 sims we want the AA geometry for no pileup sims
  //  Input::BEAM_CONFIGURATION = Input::pp_COLLISION; // for 2024 sims we want the pp geometry for no pileup sims
  //  Input::BEAM_CONFIGURATION = Input::pA_COLLISION; // for pAu sims we want the pA geometry for no pileup sims

  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  // This creates the input generator(s)
  InputInit();

  //--------------
  // Set generator specific options
  //--------------
  // can only be set after InputInit() is called

  // Simple Input generator:
  // if you run more than one of these Input::SIMPLE_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::SIMPLE)
  {
    INPUTGENERATOR::SimpleEventGenerator[0]->add_particles(simplePionPdgId, 1);
    if (Input::HEPMC || Input::EMBED)
    {
      INPUTGENERATOR::SimpleEventGenerator[0]->set_reuse_existing_vertex(true);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
    }
    else
    {
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                                PHG4SimpleEventGenerator::Uniform,
                                                                                PHG4SimpleEventGenerator::Uniform);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., simpleVertexZ);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
    }
    INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(simplePionEta, simplePionEta);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(simplePionPhi, simplePionPhi);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(simplePionPt, simplePionPt);
    std::cout << "Fun4All_G4_sPHENIX: fixed Simple pion pdg=" << simplePionPdgId
              << " pt=" << simplePionPt
              << " eta=" << simplePionEta
              << " phi=" << simplePionPhi
              << " vertex=(0,0," << simpleVertexZ << ")" << std::endl;
  }
  // Upsilons
  // if you run more than one of these Input::UPSILON_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::UPSILON)
  {
    INPUTGENERATOR::VectorMesonGenerator[0]->add_decay_particles("e", 0);
    INPUTGENERATOR::VectorMesonGenerator[0]->set_rapidity_range(-1, 1);
    INPUTGENERATOR::VectorMesonGenerator[0]->set_pt_range(0., 10.);
    // Y species - select only one, last one wins
    INPUTGENERATOR::VectorMesonGenerator[0]->set_upsilon_1s();
    if (Input::HEPMC || Input::EMBED)
    {
      INPUTGENERATOR::VectorMesonGenerator[0]->set_reuse_existing_vertex(true);
      INPUTGENERATOR::VectorMesonGenerator[0]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
    }
  }
  // particle gun
  // if you run more than one of these Input::GUN_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::GUN)
  {
    INPUTGENERATOR::Gun[0]->AddParticle("pi-", 0, 1, 0);
    INPUTGENERATOR::Gun[0]->set_vtx(0, 0, 0);
  }

  // pythia6
  if (Input::PYTHIA6)
  {
    //! Nominal collision geometry is selected by Input::BEAM_CONFIGURATION
    Input::ApplysPHENIXBeamParameter(INPUTGENERATOR::Pythia6);
  }
  // pythia8
  if (Input::PYTHIA8)
  {
    //! Nominal collision geometry is selected by Input::BEAM_CONFIGURATION
    Input::ApplysPHENIXBeamParameter(INPUTGENERATOR::Pythia8);
  }

  //--------------
  // Set Input Manager specific options
  //--------------
  // can only be set after InputInit() is called

  if (Input::HEPMC)
  {
    if (useBeamVertex)
    {
      //! Nominal collision geometry is selected by Input::BEAM_CONFIGURATION
      Input::ApplysPHENIXBeamParameter(INPUTMANAGER::HepMCInputManager);
    }
    else
    {
      std::cout << "Fun4All_G4_sPHENIX: keeping HepMC input vertices unchanged" << std::endl;
    }

    // optional overriding beam parameters
    //INPUTMANAGER::HepMCInputManager->set_vertex_distribution_width(100e-4, 100e-4, 8, 0);  //optional collision smear in space, time
    //    INPUTMANAGER::HepMCInputManager->set_vertex_distribution_mean(0,0,0,0);//optional collision central position shift in space, time
    // //optional choice of vertex distribution function in space, time
    //INPUTMANAGER::HepMCInputManager->set_vertex_distribution_function(PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus);
    //! embedding ID for the event
    //! positive ID is the embedded event of interest, e.g. jetty event from pythia
    //! negative IDs are backgrounds, .e.g out of time pile up collisions
    //! Usually, ID = 0 means the primary Au+Au collision background
    //INPUTMANAGER::HepMCInputManager->set_embedding_id(Input::EmbedID);
    if (Input::PILEUPRATE > 0)
    {
      // Copy vertex settings from foreground hepmc input
      INPUTMANAGER::HepMCPileupInputManager->CopyHelperSettings(INPUTMANAGER::HepMCInputManager);
      // and then modify the ones you want to be different
      // INPUTMANAGER::HepMCPileupInputManager->set_vertex_distribution_width(100e-4,100e-4,8,0);
    }
  }
  if (Input::PILEUPRATE > 0)
  {
    //! Nominal collision geometry is selected by Input::BEAM_CONFIGURATION
    Input::ApplysPHENIXBeamParameter(INPUTMANAGER::HepMCPileupInputManager);
  }
  // register all input generators with Fun4All
  InputRegister();

  if (! Input::READHITS)
  {
    rc->set_IntFlag("RUNNUMBER",1);

    SyncReco *sync = new SyncReco();
    se->registerSubsystem(sync);

    HeadReco *head = new HeadReco();
    se->registerSubsystem(head);
  }
// Flag Handler is always needed to read flags from input (if used)
// and update our rc flags with them. At the end it saves all flags
// again on the DST in the Flags node under the RUN node
  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  // set up production relatedstuff
  //   Enable::PRODUCTION = true;

  //======================
  // Write the DST
  //======================

  Enable::DSTOUT = false;
  Enable::DSTOUT_COMPRESS = false;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;

  auto makeOutputPath = [](const std::string &dir, const std::string &file)
  {
    if (file.empty() || (!file.empty() && file.front() == '/'))
    {
      return file;
    }
    if (dir.empty() || dir == ".")
    {
      return file;
    }
    return dir + "/" + file;
  };
  if (!outdir.empty() && outdir != ".")
  {
    gSystem->mkdir(outdir.c_str(), true);
  }
  const std::string tpcTruthTreeFile = makeOutputPath(outdir, outputFile);

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  //  Enable::DSTREADER = true;

  // turn the display on (default off)
   //Enable::DISPLAY = true;

  //======================
  // What to run
  //======================

  // QA, main switch
  Enable::QA = false;

  // Global options (enabled for all enables subsystems - if implemented)
  //  Enable::ABSORBER = true;
  //  Enable::OVERLAPCHECK = true;
  //  Enable::VERBOSITY = 1;

  // Enable::MBD = true;
  // Enable::MBD_SUPPORT = true; // save hist in MBD/BBC support structure
  // Enable::MBDRECO = Enable::MBD && true;
  Enable::MBDFAKE = false;  // Smeared vtx and t0, use if you don't want real MBD/BBC in simulation

  Enable::PIPE = true;
  Enable::PIPE_ABSORBER = false;

  // central tracking
  // ActsGeometry expects the complete tracking geometry. Build all tracker
  // volumes for TPC hitset jobs, but digitize only the TPC below.
  Enable::MVTX = !trkrHitSetDstFile.empty();
  Enable::MVTX_CELL = false;
  Enable::MVTX_CLUSTER = false;
  Enable::MVTX_QA = Enable::MVTX_CLUSTER && Enable::QA && true;

  Enable::INTT = !trkrHitSetDstFile.empty();
//  Enable::INTT_ABSORBER = true; // enables layerwise support structure readout
//  Enable::INTT_SUPPORT = true; // enable global support structure readout
  Enable::INTT_CELL = false;
  Enable::INTT_CLUSTER = false;
  Enable::INTT_QA = Enable::INTT_CLUSTER && Enable::QA && true;

  Enable::TPC = true;
  Enable::TPC_ABSORBER = false;
  // TPC_Cells creates TRKR_HITSET and TRKR_HITTRUTHASSOC.  Keep it off for
  // the original truth-point-only workflow unless a hitset DST is requested.
  Enable::TPC_CELL = !trkrHitSetDstFile.empty();
  Enable::TPC_CLUSTER = false;
  Enable::TPC_QA = Enable::TPC_CLUSTER && Enable::QA && true;

  Enable::MICROMEGAS = !trkrHitSetDstFile.empty();
  Enable::MICROMEGAS_CELL = false;
  Enable::MICROMEGAS_CLUSTER = false;
  Enable::MICROMEGAS_QA = Enable::MICROMEGAS_CLUSTER && Enable::QA && true;

  Enable::TRACKING_TRACK = false;
  Enable::GLOBAL_RECO = false;
  Enable::TRACKING_EVAL = false;
  Enable::TRACKING_QA = Enable::TRACKING_TRACK && Enable::QA && true;

  // only do track matching if TRACKINGTRACK is also used
  Enable::TRACK_MATCHING = Enable::TRACKING_TRACK && false;
  Enable::TRACK_MATCHING_TREE = Enable::TRACK_MATCHING && false;
  Enable::TRACK_MATCHING_TREE_CLUSTERS = Enable::TRACK_MATCHING_TREE && false;

  //Additional tracking tools
  //Enable::TRACKING_DIAGNOSTICS = Enable::TRACKING_TRACK && true;
  //G4TRACKING::filter_conversion_electrons = true;
  // G4TRACKING::use_alignment = true;

  // enable pp mode and set extended readout time
  // TRACKING::pp_mode = true;
  // TRACKING::pp_extended_readout_time = 20000;

  // set flags to simulate and correct TPC distortions, specify distortion and correction files
  //G4TPC::ENABLE_STATIC_DISTORTIONS = true;
  //G4TPC::static_distortion_filename = std::string("/sphenix/user/rcorliss/distortion_maps/2023.02/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX180961051_0.distortion_map.hist.root");
  //G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  //G4TPC::static_correction_filename = std::string("/sphenix/user/rcorliss/distortion_maps/2023.02/Summary_hist_mdc2_UseFieldMaps_AA_smoothed_average.correction_map.hist.root");
  //G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;

  //  cemc electronics + thin layer of W-epoxy to get albedo from cemc
  //  into the tracking, cannot run together with CEMC
  //  Enable::CEMCALBEDO = true;

  Enable::CEMC = false;
  Enable::CEMC_ABSORBER = false;
  Enable::CEMC_CELL = false;
  Enable::CEMC_TOWER = false;
  Enable::CEMC_CLUSTER = false;
  Enable::CEMC_EVAL = false;
  Enable::CEMC_QA = Enable::CEMC_CLUSTER && Enable::QA && true;

  Enable::HCALIN = false;
  Enable::HCALIN_ABSORBER = false;
  Enable::HCALIN_CELL = false;
  Enable::HCALIN_TOWER = false;
  Enable::HCALIN_CLUSTER = false;
  Enable::HCALIN_EVAL = false;
  Enable::HCALIN_QA = Enable::HCALIN_CLUSTER && Enable::QA && true;

  Enable::MAGNET = true;
  Enable::MAGNET_ABSORBER = false;

  Enable::HCALOUT = false;
  Enable::HCALOUT_ABSORBER = false;
  Enable::HCALOUT_CELL = false;
  Enable::HCALOUT_TOWER = false;
  Enable::HCALOUT_CLUSTER = false;
  Enable::HCALOUT_EVAL = false;
  Enable::HCALOUT_QA = Enable::HCALOUT_CLUSTER && Enable::QA && true;

  Enable::EPD = false;
  Enable::EPD_TILE = false;

  Enable::BEAMLINE = false;
  //  Enable::BEAMLINE_ABSORBER = true;  // makes the beam line magnets sensitive volumes
  //  Enable::BEAMLINE_BLACKHOLE = true; // turns the beamline magnets into black holes
  Enable::ZDC = false;
  //  Enable::ZDC_ABSORBER = true;
  //  Enable::ZDC_SUPPORT = true;
  Enable::ZDC_TOWER = false;
  Enable::ZDC_EVAL = false;

  //! forward flux return plug door. Out of acceptance and off by default.
  //Enable::PLUGDOOR = true;
  Enable::PLUGDOOR_ABSORBER = true;

 //Enable::GLOBAL_FASTSIM = true;

  //Enable::KFPARTICLE = true;
  //Enable::KFPARTICLE_VERBOSITY = 1;
  //Enable::KFPARTICLE_TRUTH_MATCH = true;
  //Enable::KFPARTICLE_SAVE_NTUPLE = true;

  Enable::CALOTRIGGER = Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER && false;

  Enable::JETS = false;
  Enable::JETS_EVAL = false;
  Enable::JETS_QA = Enable::JETS && Enable::QA && true;

  // HI Jet Reco for p+Au / Au+Au collisions (default is false for
  // single particle / p+p-only simulations, or for p+Au / Au+Au
  // simulations which don't particularly care about jets)
  Enable::HIJETS = Enable::JETS && Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER && false;

  // 3-D topoCluster reconstruction, potentially in all calorimeter layers
  Enable::TOPOCLUSTER = Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER && false;
  // particle flow jet reconstruction - needs topoClusters!
  Enable::PARTICLEFLOW = Enable::TOPOCLUSTER && true;
  // centrality reconstruction
  Enable::CENTRALITY = false;

  // new settings using Enable namespace in GlobalVariables.C
  Enable::BLACKHOLE = true;
  //Enable::BLACKHOLE_SAVEHITS = false; // turn off saving of bh hits
  //Enable::BLACKHOLE_FORWARD_SAVEHITS = false; // disable forward/backward hits
  //BlackHoleGeometry::visible = true;

  // run user provided code (from local G4_User.C)
  //Enable::USER = true;

  //===============
  // conditions DB flags
  //===============
  Enable::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG",CDB::global_tag);
  // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",CDB::timestamp);
  //---------------
  // World Settings
  //---------------
  //  G4WORLD::PhysicsList = "FTFP_BERT"; //FTFP_BERT_HP best for calo
  //  G4WORLD::WorldMaterial = "G4_AIR"; // set to G4_GALACTIC for material scans

  //---------------
  // Magnet Settings
  //---------------

  if (const char *uniform_field_text = gSystem->Getenv("SIM_G4_UNIFORM_BFIELD_T"))
  {
    try
    {
      const double uniform_field_t = std::stod(uniform_field_text);
      if (!std::isfinite(uniform_field_t) || uniform_field_t == 0.0)
      {
        throw std::out_of_range("field must be finite and nonzero");
      }
      G4MAGNET::magfield = uniform_field_text;
      std::cout << "Fun4All_G4_sPHENIX: fixed uniform magnetic field="
                << uniform_field_t << " T" << std::endl;
    }
    catch (const std::exception &error)
    {
      std::cout << "Fun4All_G4_sPHENIX: invalid SIM_G4_UNIFORM_BFIELD_T='"
                << uniform_field_text << "': " << error.what() << std::endl;
      return 1;
    }
  }

  //  G4MAGNET::magfield =  std::string(getenv("CALIBRATIONROOT"))+ std::string("/Field/Map/sphenix3dbigmapxyz.root");  // default map from the calibration database
  //  G4MAGNET::magfield = "1.5"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
//  G4MAGNET::magfield_rescale = 1.;  // make consistent with expected Babar field strength of 1.4T

  //---------------
  // Pythia Decayer
  //---------------
  // list of decay types in
  // $OFFLINE_MAIN/include/g4decayer/EDecayType.hh
  // default is All:
  // G4P6DECAYER::decayType = EDecayType::kAll;

  // Initialize the selected subsystems
  G4Init();
  // Keep the standard sPHENIX tracking layer numbering even though this
  // truth-point job does not build the MVTX/INTT detector volumes.
  G4MVTX::n_maps_layer = 3;
  G4INTT::n_intt_layer = 4;

  //---------------------
  // GEANT4 Detector description
  //---------------------
  if (!Input::READHITS)
  {
    G4Setup();
  }

  //------------------
  // Detector Division
  //------------------

  if ((Enable::MBD && Enable::MBDRECO) || Enable::MBDFAKE) Mbd_Reco();

  if (Enable::CEMC_CELL) CEMC_Cells();

  if (Enable::HCALIN_CELL) HCALInner_Cells();

  if (Enable::HCALOUT_CELL) HCALOuter_Cells();

  if (Enable::TPC_CELL)
  {
    TPC_Cells();
    // PHG4TpcElectronDrift converts global positions through ActsGeometry.
    // The truth-point-only path does not otherwise need to build this node.
    ACTSGEOM::ActsGeomInit();
  }

  // An empty outputFile is useful for jobs that only want TRKR_HITSET.
  if (!outputFile.empty() || !truthPointDstFile.empty())
  {
    auto *tpcTruthPoints = new PHG4TpcTruthPointBuilder();
    tpcTruthPoints->set_output_node("G4HIT_TPC_TRUECLUSTER");
    // Keep every crossing for low-pT loopers that revisit the same layer.
    tpcTruthPoints->set_max_intersections_per_track_layer(0);
    tpcTruthPoints->Verbosity(1);
    se->registerSubsystem(tpcTruthPoints);

    if (!outputFile.empty())
    {
      auto *tpcTruthTree = new PHG4TpcTruthPointTree("PHG4TpcTruthPointTree", tpcTruthTreeFile);
      tpcTruthTree->set_truth_point_node("G4HIT_TPC_TRUECLUSTER");
      tpcTruthTree->Verbosity(1);
      se->registerSubsystem(tpcTruthTree);
    }
  }

  //-----------------------------
  // CEMC towering and clustering
  //-----------------------------

  if (Enable::CEMC_TOWER) CEMC_Towers();
  if (Enable::CEMC_CLUSTER) CEMC_Clusters();

  //--------------
  // EPD tile reconstruction
  //--------------

  if (Enable::EPD_TILE) EPD_Tiles();

  //-----------------------------
  // HCAL towering and clustering
  //-----------------------------

  if (Enable::HCALIN_TOWER) HCALInner_Towers();
  if (Enable::HCALIN_CLUSTER) HCALInner_Clusters();

  if (Enable::HCALOUT_TOWER) HCALOuter_Towers();
  if (Enable::HCALOUT_CLUSTER) HCALOuter_Clusters();

  // if enabled, do topoClustering early, upstream of any possible jet reconstruction
  if (Enable::TOPOCLUSTER) TopoClusterReco();

  // This stripped-down macro writes ideal TPC truth points directly from
  // G4HIT_TPC, so it intentionally skips cellization, clustering, and track
  // reconstruction.


  //-----------------
  // Global Vertexing
  //-----------------

  if (Enable::GLOBAL_RECO && Enable::GLOBAL_FASTSIM)
  {
    std::cout << "You can only enable Enable::GLOBAL_RECO or Enable::GLOBAL_FASTSIM, not both" << std::endl;
    gSystem->Exit(1);
  }
  if (Enable::GLOBAL_RECO)
  {
    Global_Reco();
  }
  else if (Enable::GLOBAL_FASTSIM)
  {
    Global_FastSim();
  }

  //-----------------
  // Centrality Determination
  //-----------------

  if (Enable::CENTRALITY)
  {
      Centrality();
  }

  //-----------------
  // Calo Trigger Simulation
  //-----------------

  if (Enable::CALOTRIGGER)
  {
    CaloTrigger_Sim();
  }

  //---------
  // Jet reco
  //---------

  if (Enable::JETS) Jet_Reco();
  if (Enable::HIJETS) HIJetReco();

  if (Enable::PARTICLEFLOW) ParticleFlow();

  //----------------------
  // Simulation evaluation
  //----------------------
  std::string outputroot = outputFile;
  std::string remove_this = ".root";
  size_t pos = outputroot.find(remove_this);
  if (pos != std::string::npos)
  {
    outputroot.erase(pos, remove_this.length());
  }

  if (Enable::CEMC_EVAL) CEMC_Eval(outputroot + "_g4cemc_eval.root");

  if (Enable::HCALIN_EVAL) HCALInner_Eval(outputroot + "_g4hcalin_eval.root");

  if (Enable::HCALOUT_EVAL) HCALOuter_Eval(outputroot + "_g4hcalout_eval.root");

  if (Enable::JETS_EVAL) Jet_Eval(outputroot + "_g4jet_eval.root");

  if (Enable::DSTREADER) G4DSTreader(outputroot + "_DSTReader.root");



  if (Enable::USER) UserAnalysisInit();

  //======================
  // Run KFParticle on evt
  //======================
  if (Enable::KFPARTICLE && Input::UPSILON) KFParticle_Upsilon_Reco();
  if (Enable::KFPARTICLE && Input::DZERO) KFParticle_D0_Reco();

  //----------------------
  // Standard QAs
  //----------------------

  if (Enable::CEMC_QA) CEMC_QA();
  if (Enable::HCALIN_QA) HCALInner_QA();
  if (Enable::HCALOUT_QA) HCALOuter_QA();

  if (Enable::JETS_QA) Jet_QA();

  if (Enable::TRACKING_QA && Enable::CEMC_QA && Enable::HCALIN_QA && Enable::HCALOUT_QA) QA_G4CaloTracking();

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  if (Enable::PRODUCTION)
  {
    Production_CreateOutputDir();
  }

  if (Enable::DSTOUT)
  {
    std::string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    if (Enable::DSTOUT_COMPRESS)
    {
      ShowerCompress();
      DstCompress(out);
    }
    se->registerOutputManager(out);
  }

  if (!truthPointDstFile.empty())
  {
    auto *truthPointOut = new Fun4AllDstOutputManager("TPC_TRUTH_POINT_DST", truthPointDstFile);
    truthPointOut->AddNode("EventHeader");
    truthPointOut->AddNode("G4TruthInfo");
    truthPointOut->AddNode("G4HIT_TPC_TRUECLUSTER");
    se->registerOutputManager(truthPointOut);
  }

  if (!trkrHitSetDstFile.empty())
  {
    auto *trkrHitSetOut = new Fun4AllDstOutputManager("TPC_TRKR_HITSET_DST", trkrHitSetDstFile);
    trkrHitSetOut->AddNode("TRKR_HITSET");
    trkrHitSetOut->AddNode("TRKR_HITTRUTHASSOC");
    se->registerOutputManager(trkrHitSetOut);
  }
  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY)
  {
    DisplayOn();

    gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
    gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "You are in event display mode. Run one event with" << std::endl;
    std::cout << "se->run(1)" << std::endl;
    std::cout << "Run Geant4 command with following examples" << std::endl;
    gROOT->ProcessLine("displaycmd()");

    return 0;
  }

  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0)
  {
    return 0;
  }
  // if we run the particle generator and use 0 it'll run forever
  // for embedding it runs forever if the repeat flag is set
  if (nEvents == 0 && !Input::HEPMC && !Input::READHITS && INPUTEMBED::REPEAT)
  {
    std::cout << "using 0 for number of events is a bad idea when using particle generators" << std::endl;
    std::cout << "it will run forever, so I just return without running anything" << std::endl;
    return 0;
  }

  se->skip(skip);
  se->run(nEvents);
  //  se->PrintTimer();

  //-----
  // QA output
  //-----

  if (Enable::QA) QA_Output(outputroot + "_qa.root");

  //-----
  // Exit
  //-----

  CDBInterface::instance()->Print(); // print used DB files
  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  if (Enable::PRODUCTION)
  {
    Production_MoveOutput();
  }

  gSystem->Exit(0);
  return 0;
}
#endif
