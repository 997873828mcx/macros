// -----------------------------------------------------------------------------
// Run the lightweight KshortReconstruction AP-tree module on reconstructed
// real-data track DSTs.
//
// Example:
//   root -l -b -q 'Fun4All_KshortAP.C(1000, "dst_tracks.list", "ap_real.root")'
//
// If the input is a single DST instead of a list:
//   root -l -b -q 'Fun4All_KshortAP.C(1000, "DST_TRKR_TRACKS.root", "ap_real.root", false)'
// -----------------------------------------------------------------------------

#include <GlobalVariables.C>

#include <G4_ActsGeom.C>
#include <G4_Magnet.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/recoConsts.h>

#include <trackingdiagnostics/KshortReconstruction.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)

namespace
{
std::string first_input_file(const std::string& input, const bool input_is_list)
{
  if (!input_is_list)
  {
    return input;
  }

  std::ifstream in(input);
  std::string line;
  while (std::getline(in, line))
  {
    if (!line.empty() && line[0] != '#')
    {
      return line;
    }
  }
  return "";
}
}  // namespace

void Fun4All_KshortAP(const int nEvents = 1000,
                      const std::string& input = "dst_tracks.list",
                      const std::string& output = "ap_realdata.root",
                      const bool inputIsList = true,
                      const int nSkip = 0,
                      const bool requireMvtx = false,
                      const double trackPtCut = 0.05,
                      const int runNumberOverride = 0)
{
  auto* se = Fun4AllServer::instance();
  se->Verbosity(1);

  std::string first_file = first_input_file(input, inputIsList);
  if (first_file.empty())
  {
    std::cerr << "Fun4All_KshortAP: no input file found in " << input << std::endl;
    gSystem->Exit(1);
  }

  int runnumber = runNumberOverride;
  if (runnumber <= 0)
  {
    runnumber = Fun4AllUtils::GetRunSegment(first_file).first;
  }
  if (runnumber <= 0)
  {
    std::cerr << "Fun4All_KshortAP: could not infer run number from " << first_file
              << ". Pass runNumberOverride as the 8th argument." << std::endl;
    gSystem->Exit(1);
  }

  auto* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");

  TpcReadoutInit(runnumber);
  std::cout << "Fun4All_KshortAP: run " << runnumber
            << " samples " << TRACKING::reco_tpc_maxtime_sample
            << " pre " << TRACKING::reco_tpc_time_presample
            << " vdrift " << G4TPC::tpc_drift_velocity_reco
            << std::endl;

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  auto* ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TrackingInit();
  TRACKING::pp_mode = true;

  auto* trackin = new Fun4AllDstInputManager("TrackInput");
  if (inputIsList)
  {
    trackin->AddListFile(input);
  }
  else
  {
    trackin->fileopen(input);
  }
  se->registerInputManager(trackin);

  auto* kshort = new KshortReconstruction("KshortReconstruction_AP");
  kshort->set_output_file(output);
  kshort->setRequireMVTX(requireMvtx);
  kshort->setTrackPtCut(trackPtCut);
  kshort->setPairDCACut(0.3);
  kshort->setTrackDCACut(0.01);
  kshort->setApplyInvariantPtCut(false);
  kshort->Verbosity(0);
  se->registerSubsystem(kshort);

  if (nSkip > 0)
  {
    se->skip(nSkip);
  }
  se->run(nEvents);
  se->End();
  se->PrintTimer();

  delete se;
  std::cout << "Fun4All_KshortAP: wrote " << output << std::endl;
  gSystem->Exit(0);
}
