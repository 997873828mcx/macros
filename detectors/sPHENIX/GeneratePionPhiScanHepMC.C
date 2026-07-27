#ifndef MACRO_GENERATEPIONPHISCANHEPMC_C
#define MACRO_GENERATEPIONPHISCANHEPMC_C

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/Units.h>
#pragma GCC diagnostic pop

#include <Rtypes.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(libHepMC.so)

// Write identical events containing an evenly spaced phi scan of stable pions.
// At eta=0, the requested total momentum p is also pT.
int GeneratePionPhiScanHepMC(
    const std::string &outputFile,
    const int pionPdgId = 211,
    const int nPions = 24,
    const double minMomentum = 0.2,
    const double maxMomentum = 1.0,
    const double eta = 0.0,
    const int nEvents = 1)
{
  if ((pionPdgId != 211 && pionPdgId != -211) || nPions <= 0 ||
      minMomentum <= 0 || maxMomentum < minMomentum || nEvents <= 0)
  {
    std::cout << "GeneratePionPhiScanHepMC: invalid arguments" << std::endl;
    return 1;
  }

  const auto slash = outputFile.find_last_of('/');
  if (slash != std::string::npos)
  {
    const std::string outputDir = outputFile.substr(0, slash);
    if (!outputDir.empty()) gSystem->mkdir(outputDir.c_str(), true);
  }

  constexpr double pionMass = 0.13957039;  // GeV
  constexpr double twoPi = 2.0 * M_PI;
  HepMC::IO_GenEvent writer(outputFile, std::ios::out);

  for (int eventNumber = 0; eventNumber < nEvents; ++eventNumber)
  {
    auto *event = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
    event->set_event_number(eventNumber);

    auto *vertex = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0., 0.));
    event->add_vertex(vertex);

    for (int i = 0; i < nPions; ++i)
    {
      const double fraction = (nPions == 1) ? 0.0 : static_cast<double>(i) / (nPions - 1);
      const double momentum = minMomentum + fraction * (maxMomentum - minMomentum);
      const double phi = twoPi * static_cast<double>(i) / nPions;
      const double pt = momentum / std::cosh(eta);
      const double px = pt * std::cos(phi);
      const double py = pt * std::sin(phi);
      const double pz = pt * std::sinh(eta);
      const double energy = std::sqrt(momentum * momentum + pionMass * pionMass);

      auto *pion = new HepMC::GenParticle(
          HepMC::FourVector(px, py, pz, energy), pionPdgId, 1);
      vertex->add_particle_out(pion);
    }

    writer << event;
    delete event;
  }

  std::cout << "GeneratePionPhiScanHepMC: wrote " << nEvents << " event(s), each with "
            << nPions << (pionPdgId > 0 ? " pi+" : " pi-") << ", to " << outputFile
            << std::endl;
  return 0;
}

#endif
