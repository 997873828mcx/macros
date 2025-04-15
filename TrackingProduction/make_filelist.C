#ifndef make_filelist_C
#define make_filelist_C

//____________________________________________________________________
std::vector<std::string> make_filelist( const std::string& tag = "ana464_nocdbtag_v001", const int runnumber = 53756, const int segment = 0)
{

  const std::vector<std::string> detector_tags = {
    "DST_STREAMING_EVENT_INTT0",
    "DST_STREAMING_EVENT_INTT1",
    "DST_STREAMING_EVENT_INTT2",
    "DST_STREAMING_EVENT_INTT3",
    "DST_STREAMING_EVENT_INTT4",
    "DST_STREAMING_EVENT_INTT5",
    "DST_STREAMING_EVENT_INTT6",
    "DST_STREAMING_EVENT_INTT7",
    "DST_STREAMING_EVENT_MVTX0",
    "DST_STREAMING_EVENT_MVTX1",
    "DST_STREAMING_EVENT_MVTX2",
    "DST_STREAMING_EVENT_MVTX3",
    "DST_STREAMING_EVENT_MVTX4",
    "DST_STREAMING_EVENT_MVTX5",
    "DST_STREAMING_EVENT_TPC00",
    "DST_STREAMING_EVENT_TPC01",
    "DST_STREAMING_EVENT_TPC02",
    "DST_STREAMING_EVENT_TPC03",
    "DST_STREAMING_EVENT_TPC04",
    "DST_STREAMING_EVENT_TPC05",
    "DST_STREAMING_EVENT_TPC06",
    "DST_STREAMING_EVENT_TPC07",
    "DST_STREAMING_EVENT_TPC08",
    "DST_STREAMING_EVENT_TPC09",
    "DST_STREAMING_EVENT_TPC10",
    "DST_STREAMING_EVENT_TPC11",
    "DST_STREAMING_EVENT_TPC12",
    "DST_STREAMING_EVENT_TPC13",
    "DST_STREAMING_EVENT_TPC14",
    "DST_STREAMING_EVENT_TPC15",
    "DST_STREAMING_EVENT_TPC16",
    "DST_STREAMING_EVENT_TPC17",
    "DST_STREAMING_EVENT_TPC18",
    "DST_STREAMING_EVENT_TPC19",
    "DST_STREAMING_EVENT_TPC20",
    "DST_STREAMING_EVENT_TPC21",
    "DST_STREAMING_EVENT_TPC22",
    "DST_STREAMING_EVENT_TPC23",
    "DST_STREAMING_EVENT_TPOT"
  };

  // get run range
  const int range_begin = int( runnumber/100 ) * 100;
  const int range_end = int( runnumber/100 + 1 ) * 100;

  std::vector<std::string> filelist;

  // loop over detector tags
  for( const auto& detector_tag:detector_tags )
  {
    std::string filename = Form( "/sphenix/lustre01/sphnxpro/production/run2pp/physics/%s/%s/run_%08i_%08i/dst/%s_run2pp_%s-%08i-%05i.root",
      tag.c_str(), detector_tag.c_str(),
      range_begin, range_end,
      detector_tag.c_str(), tag.c_str(),
      runnumber, segment );

    if( ifstream(filename.c_str()).good() )
    {
      std::cout << "make_filelist - adding: " << filename << std::endl;
      filelist.emplace_back(filename);
    }
  }

  if( filelist.empty() )
  { std::cout << "make_filelist - no files found" << std::endl; }

  return filelist;
}

#endif
