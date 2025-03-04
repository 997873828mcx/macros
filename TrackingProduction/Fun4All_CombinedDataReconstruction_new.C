#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

// Micromegas
#include <micromegas/MicromegasCombinedDataDecoder.h>
#include <micromegas/MicromegasCombinedDataEvaluation.h>
#include <micromegas/MicromegasClusterizer.h>

// INTT
#include <intt/InttCombinedRawDataDecoder.h>
#include <intt/InttClusterizer.h>

// TPC
#include <tpc/TpcCombinedRawDataUnpacker.h>

#include <tpccalib/PHTpcResiduals.h>

#include <trackingdiagnostics/TrackResiduals.h>

// local macros
#include "G4Setup_sPHENIX.C"
#include "G4_Global.C"

#include "Trkr_RecoInit.C"
#include "Trkr_TpcReadoutInit.C"

#include "Trkr_Clustering.C"
#include "Trkr_Reco.C"

R__LOAD_LIBRARY(libfun4all.so)

R__LOAD_LIBRARY(libmicromegas.so)

// tag = "ana464_nocdbtag_v001/"

#include "make_filelist.C"

int Fun4All_CombinedDataReconstruction_new(
    const int nEvents = 1,
    const int nSkipEvents = 0,
    const char *tag = "ana464_nocdbtag_v001",
    const int runnumber = 53534,
    const int segment = 0,
    const char *outputFile = "dst_eval-00053534-0000_corrected.root",
    const char *residualFile = "TrackResiduals-00053534-0000_corrected.root")
{
    // print inputs
    std::cout << "Fun4All_CombinedDataReconstruction - nEvents: " << nEvents << std::endl;
    std::cout << "Fun4All_CombinedDataReconstruction - nSkipEvents: " << nSkipEvents << std::endl;
    std::cout << "Fun4All_CombinedDataReconstruction - tag: " << tag << std::endl;
    std::cout << "Fun4All_CombinedDataReconstruction - runnumber: " << runnumber << std::endl;
    std::cout << "Fun4All_CombinedDataReconstruction - segment: " << segment << std::endl;
    std::cout << "Fun4All_CombinedDataReconstruction - outputFile: " << outputFile << std::endl;
    std::cout << "Fun4All_CombinedDataReconstruction - residualFile: " << residualFile << std::endl;

    std::string outputFileName = "outputFile_" + std::to_string(runnumber) + "_" + std::to_string(segment);
    TRACKING::pp_mode = true;
    G4TRACKING::SC_CALIBMODE = false;
    G4TRACKING::convert_seeds_to_svtxtracks = false;

    // condition database
    Enable::CDB = true;

    // reco const
    auto rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", runnumber);
    rc->set_IntFlag("RUNSEGMENT", segment);
    rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
    rc->set_uint64Flag("TIMESTAMP", runnumber);

    // try use hard coded geometry file
    // /cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/Tracking_Geometry/./87/6a/876abfa36f58e100f29e6641f9ff4f5b_tracking_geometry_50nsclock.root
    TpcReadoutInit(runnumber);

    // TPC unpacker t0
    int tpcunpacker_t0 = 0;

    if (false)
    {

        // need to apply proper scaling to match new TPC FEE clock frequency
        G4TPC::tpc_drift_velocity_reco = G4TPC::tpc_drift_velocity_reco * 53.0 / 50.0373;

        // reassign to ZCrossing correction
        TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;

        // unpacker t0, from bade
        tpcunpacker_t0 = -4;
    }
    else
    {

        // values from Tony
        G4TPC::tpc_drift_velocity_reco = 0.00725;

        // reassign to ZCrossing correction
        TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;

        // values from Tony
        tpcunpacker_t0 = -9;
    }

    // printout
    std::cout << "Fun4All_CombinedDataReconstruction - samples: " << TRACKING::reco_tpc_maxtime_sample << std::endl;
    std::cout << "Fun4All_CombinedDataReconstruction - pre: " << TRACKING::reco_tpc_time_presample << std::endl;
    std::cout << "Fun4All_CombinedDataReconstruction - vdrift: " << G4TPC::tpc_drift_velocity_reco << std::endl;

    ACTSGEOM::mvtxMisalignment = 100;
    ACTSGEOM::inttMisalignment = 100;
    ACTSGEOM::tpotMisalignment = 100;

    // server
    auto se = Fun4AllServer::instance();
    se->Verbosity(0);

    PHRandomSeed::Verbosity(1);

    // tracking geometry
    std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
    std::cout << "Fun4All_CombinedDataReconstruction_hp - geofile: " << geofile << std::endl;

    auto ingeo = new Fun4AllRunNodeInputManager("GeoIn");
    ingeo->AddFile(geofile);
    se->registerInputManager(ingeo);

    // distortion correction
    if (false)
    {
        // module edge corrections
        G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;

        // static distortions
        G4TPC::ENABLE_STATIC_CORRECTIONS = true;

        // average distortions
        G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;
    }

    std::string edgeStatus = (G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS ? "edgeOn" : "edgeOff");
    std::string staticStatus = (G4TPC::ENABLE_STATIC_CORRECTIONS ? "staticOn" : "staticOff");
    std::string seedType = (G4TRACKING::convert_seeds_to_svtxtracks ? "seeds" : "acts");

    std::time_t now = std::time(nullptr);
    std::tm tm = *std::localtime(&now);
    char dateStr[9];
    std::strftime(dateStr, sizeof(dateStr), "%m%d", &tm);

    // tpc zero suppression
    TRACKING::tpc_zero_supp = true;

    G4MAGNET::magfield_rescale = 1;
    TrackingInit();

    // input managers
    {
        const auto filelist = make_filelist(tag, runnumber, segment);
        for (size_t i = 0; i < filelist.size(); ++i)
        {
            const auto &inputfile = filelist[i];
            auto in = new Fun4AllDstInputManager(Form("DSTin_%zu", i));
            in->fileopen(inputfile);
            se->registerInputManager(in);
        }
    }

    // hit unpackers
    for (int felix = 0; felix < 6; felix++)
    {
        Mvtx_HitUnpacking(std::to_string(felix));
    }

    for (int server = 0; server < 8; server++)
    {
        Intt_HitUnpacking(std::to_string(server));
    }

    {
        // TPC unpacking
        for (int ebdc = 0; ebdc < 24; ebdc++)
        {
            const std::string ebdc_string = (Form("%02i", ebdc));

            auto tpcunpacker = new TpcCombinedRawDataUnpacker("TpcCombinedRawDataUnpacker" + ebdc_string);
            tpcunpacker->useRawHitNodeName("TPCRAWHIT_" + ebdc_string);
            tpcunpacker->set_presampleShift(TRACKING::reco_tpc_time_presample);
            tpcunpacker->set_t0(tpcunpacker_t0);
            if (TRACKING::tpc_zero_supp)
            {
                tpcunpacker->ReadZeroSuppressedData();
            }

            tpcunpacker->doBaselineCorr(TRACKING::tpc_baseline_corr);
            se->registerSubsystem(tpcunpacker);
        }
    }

    {
        // micromegas unpacking
        auto tpotunpacker = new MicromegasCombinedDataDecoder;
        const auto calibrationFile = CDBInterface::instance()->getUrl("TPOT_Pedestal");
        tpotunpacker->set_calibration_file(calibrationFile);
        tpotunpacker->set_sample_max(1024);
        // tpotunpacker->set_hot_channel_map_file("Calibrations/TPOT_HotChannels-00041989-0000.root" );
        se->registerSubsystem(tpotunpacker);
    }

    Mvtx_Clustering();
    Intt_Clustering();

    // tpc clustering
    {

        TString outfileTpcClusterizer = TString::Format("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/%s_clusterizer_%s_%s_%s_%s.root",
                                                        outputFileName.c_str(),
                                                        edgeStatus.c_str(),
                                                        staticStatus.c_str(),
                                                        seedType.c_str(),
                                                        dateStr);
        std::string outTpcClusterizerString(outfileTpcClusterizer.Data());
        auto tpcclusterizer = new TpcClusterizer;
        tpcclusterizer->Verbosity(0);
        tpcclusterizer->set_outfileName(outTpcClusterizerString);
        tpcclusterizer->set_rawdata_reco();
        tpcclusterizer->set_sampa_tbias(0);
        se->registerSubsystem(tpcclusterizer);
    }

    Tpc_LaserEventIdentifying();

    Micromegas_Clustering();

    {
        // silicon seeding

        auto silicon_Seeding = new PHActsSiliconSeeding;
        silicon_Seeding->Verbosity(0);
        silicon_Seeding->setinttRPhiSearchWindow(1.0);
        silicon_Seeding->setinttZSearchWindow(7.0);
        silicon_Seeding->seedAnalysis(false);
        se->registerSubsystem(silicon_Seeding);
    }

    {
        auto merger = new PHSiliconSeedMerger;
        merger->Verbosity(0);
        se->registerSubsystem(merger);
    }

    double fieldstrength = std::numeric_limits<double>::quiet_NaN();
    const bool ConstField = isConstantField(G4MAGNET::magfield_tracking, fieldstrength);

    {
        // TPC seeding
        TString outfileSeed = TString::Format("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/%s_clusters_%s_%s_%s_%s.root",
                                              outputFileName.c_str(),
                                              edgeStatus.c_str(),
                                              staticStatus.c_str(),
                                              seedType.c_str(),
                                              dateStr);
        std::string outstring(outfileSeed.Data());
        auto seeder = new PHCASeeding("PHCASeeding");
        if (ConstField)
        {
            seeder->useConstBField(true);
            seeder->constBField(fieldstrength);
        }
        else
        {
            seeder->set_field_dir(-1 * G4MAGNET::magfield_rescale);
            seeder->useConstBField(false);
            seeder->magFieldFile(G4MAGNET::magfield_tracking); // to get charge sign right
        }
        seeder->Verbosity(0);
        seeder->SetLayerRange(7, 55);
        seeder->SetSearchWindow(2., 0.05);          // z-width and phi-width, default in macro at 1.5 and 0.05
        seeder->SetClusAdd_delta_window(3.0, 0.06); //  (0.5, 0.005) are default; sdzdr_cutoff, d2/dr2(phi)_cutoff
        seeder->SetMinHitsPerCluster(0);
        seeder->SetMinClustersPerTrack(3);
        seeder->useFixedClusterError(true);
        seeder->set_pp_mode(true);
        seeder->set_outfileName(outstring);
        // seeder->set_pp_mode(TRACKING::pp_mode);
        se->registerSubsystem(seeder);
    }

    {
        // expand stubs in the TPC using simple kalman filter
        auto cprop = new PHSimpleKFProp("PHSimpleKFProp");
        cprop->set_field_dir(G4MAGNET::magfield_rescale);
        if (ConstField)
        {
            cprop->useConstBField(true);
            cprop->setConstBField(fieldstrength);
        }
        else
        {
            cprop->magFieldFile(G4MAGNET::magfield_tracking);
            cprop->set_field_dir(-1 * G4MAGNET::magfield_rescale);
        }
        cprop->useFixedClusterError(true);
        cprop->set_max_window(5.);
        cprop->Verbosity(0);
        // cprop->set_pp_mode(TRACKING::pp_mode);
        cprop->set_pp_mode(true);
        se->registerSubsystem(cprop);

        if (TRACKING::pp_mode)
        {
            // for pp mode, apply preliminary distortion corrections to TPC clusters before crossing is known
            // and refit the trackseeds. Replace KFProp fits with the new fit parameters in the TPC seeds.
            auto prelim_distcorr = new PrelimDistortionCorrection;
            prelim_distcorr->set_pp_mode(TRACKING::pp_mode);
            prelim_distcorr->Verbosity(0);
            se->registerSubsystem(prelim_distcorr);
        }
    }

    {
        // matching to silicons
        auto silicon_match = new PHSiliconTpcTrackMatching;
        silicon_match->Verbosity(0);

        // wide matching windows
        //     silicon_match->set_x_search_window(2.);
        //     silicon_match->set_y_search_window(2.);
        //     silicon_match->set_z_search_window(5.);
        //     silicon_match->set_phi_search_window(0.2);
        //     silicon_match->set_eta_search_window(0.1);

        // narrow matching windows
        silicon_match->set_x_search_window(0.36);
        silicon_match->set_y_search_window(0.36);
        silicon_match->set_z_search_window(2.5);
        silicon_match->set_phi_search_window(0.014);
        silicon_match->set_eta_search_window(0.0091);
        silicon_match->set_test_windows_printout(false);

        silicon_match->set_pp_mode(TRACKING::pp_mode);
        se->registerSubsystem(silicon_match);
    }

    if (true)
    {
        // matching with micromegas
        auto mm_match = new PHMicromegasTpcTrackMatching;
        mm_match->Verbosity(0);

        mm_match->set_rphi_search_window_lyr1(3.0);
        mm_match->set_rphi_search_window_lyr2(15.0);

        mm_match->set_z_search_window_lyr1(30.0);
        mm_match->set_z_search_window_lyr2(3.0);

        mm_match->set_min_tpc_layer(38);            // layer in TPC to start projection fit
        mm_match->set_test_windows_printout(false); // used for tuning search windows only
        se->registerSubsystem(mm_match);
    }

    {
        // track fit
        se->registerSubsystem(new PHTpcDeltaZCorrection);

        // perform final track fit with ACTS
        auto actsFit = new PHActsTrkFitter;
        actsFit->Verbosity(0);
        actsFit->commissioning(G4TRACKING::use_alignment);

        // fit with Micromegas and Silicon ONLY
        actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
        actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);

        actsFit->set_pp_mode(TRACKING::pp_mode);
        actsFit->set_use_clustermover(true);
        actsFit->useActsEvaluator(false);
        actsFit->useOutlierFinder(false);
        actsFit->setFieldMap(G4MAGNET::magfield_tracking);
        se->registerSubsystem(actsFit);
    }

    {
        auto cleaner = new PHTrackCleaner();
        cleaner->Verbosity(0);
        se->registerSubsystem(cleaner);
    }

    {
        // vertex finding
        auto finder = new PHSimpleVertexFinder;
        finder->Verbosity(0);
        finder->setDcaCut(0.5);
        finder->setTrackPtCut(-99999.);
        finder->setBeamLineCut(1);
        finder->setTrackQualityCut(100);
        finder->setNmvtxRequired(3);
        finder->setOutlierPairCut(1);

        se->registerSubsystem(finder);
    }

    if (G4TRACKING::SC_CALIBMODE)
    {
        /*
         * in calibration mode, calculate residuals between TPC and fitted tracks,
         * store in dedicated structure for distortion correction
         */
        auto residuals = new PHTpcResiduals;
        residuals->setOutputfile("PhTpcResiduals.root");
        residuals->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);

        // matches Tony's analysis
        residuals->setMinPt(0.2);

        // reconstructed distortion grid size (phi, r, z)
        residuals->setGridDimensions(36, 48, 80);
        se->registerSubsystem(residuals);
    }

    if (true)
    {

        TString residoutfile = TString::Format("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/%s_resid_%s_%s_%s_%s.root",
                                               outputFileName.c_str(),
                                               edgeStatus.c_str(),
                                               staticStatus.c_str(),
                                               seedType.c_str(),
                                               dateStr);
        std::string residstring(residoutfile.Data());
        auto resid = new TrackResiduals("TrackResiduals");
        resid->outfileName(residstring);
        resid->Verbosity(0);
        resid->alignment(false);

        // resid->set_rejectLaserEvent(false);

        // adjust track map made
        if (G4TRACKING::SC_CALIBMODE)
        {
            resid->trackmapName("SvtxSiliconMMTrackMap");
        }

        resid->clusterTree();
        resid->hitTree();
        resid->convertSeeds(G4TRACKING::convert_seeds_to_svtxtracks);
        // resid->linefitAll();
        se->registerSubsystem(resid);
    }

    if (true)
    {
        auto out = new Fun4AllDstOutputManager("DSTOUT", outputFile);
        out->AddNode("TrackingEvaluator_hp::Container");
        out->AddNode("MicromegasTrackEvaluator_hp::Container");
        se->registerOutputManager(out);
    }

    // process events
    if (nSkipEvents > 0)
    {
        se->skip(nSkipEvents);
    }

    se->run(nEvents);

    // terminate
    se->End();
    se->PrintTimer();

    std::cout << "All done" << std::endl;
    delete se;
    gSystem->Exit(0);
    return 0;
}

// This function is only used to test if we can load this as root6 macro
// without running into unresolved libraries and include files
void RunFFALoadTest() {}
