#include "nain4.hh"
#include "g4-mandatory.hh"

#include "geometries/imas.hh"
#include "geometries/nema.hh"
#include "geometries/samples.hh"
#include "geometries/sipm.hh"

#include <G4RunManager.hh>
#include <G4RunManagerFactory.hh>
#include <G4SystemOfUnits.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4VisManager.hh>

#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4OpticalPhysics.hh>
#include <Randomize.hh>

#include <memory>
#include <string>

using std::make_unique;
using std::unique_ptr;

int main(int argc, char** argv) {
  // Detect interactive mode (if no arguments) and define UI session
  auto ui = argc == 1
    ? make_unique<G4UIExecutive>(argc, argv)
    : unique_ptr <G4UIExecutive>{nullptr};

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);

  // Construct the default run manager
  auto run_manager = unique_ptr<G4RunManager>
    {G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial)};

  // For use with phantom_in_cylinder
  auto phantom = a_nema_phantom();

  int last_interesting_event = 0;
  size_t count_interesting_event = 0;
  size_t count_511s              = 0;
  size_t count_gammas            = 0;
  size_t events_with_electrons = 0;  bool got_electron = false;
  size_t events_with_optphots  = 0;  bool got_optphot  = false;
  // clang-format off
  n4::sensitive_detector::process_hits_fn make_noise = [&last_interesting_event, &count_interesting_event,
                                                        &events_with_electrons, &got_electron,
                                                        &events_with_optphots , &got_optphot,
                                                        &count_511s, &count_gammas](G4Step* step) {
    auto pos      = step -> GetPostStepPoint() -> GetPosition();
    auto track    = step -> GetTrack();
    auto energy   = track -> GetKineticEnergy();
    auto particle = track -> GetParticleDefinition();
    auto name     = particle -> GetParticleName();
    auto id       = track -> GetTrackID();
    auto event_id = n4::event_number();
    auto tg = track -> GetGlobalTime();

    auto vol_pre_name = step -> GetPreStepPoint()  -> GetTouchable() -> GetVolume() -> GetName();
    auto vol_postname = step -> GetPostStepPoint() -> GetTouchable() -> GetVolume() -> GetName();

    // proc_name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

    if (event_id > last_interesting_event) { count_interesting_event++; }
    last_interesting_event = event_id;

    if (name == "gamma") {
      count_gammas++;
      if (energy > 510.99 * keV && energy < 511.01 * keV) {
        count_511s++;
      }
    }

    if (name == "e-") {
      if (!got_electron) { events_with_electrons++; }
      got_electron = true;
    }

    if (name == "opticalphoton") {
      if (!got_optphot) { events_with_optphots++; }
      got_electron = true;
    }
    std::cout << std::setw (4) << events_with_electrons << " /"
              << std::setw (4) << events_with_optphots << " /"
              << std::setw (4) << round(100 - (100.0 * count_511s / count_gammas) ) << "% /"
              << std::setw (4) << count_interesting_event << " /"
              << std::setw (4) << event_id << ' '
              << std::setw(15) << name << ' '
              << std::setw (4) << id << ' '
              << std::setw (4) << round(energy / keV) << " keV " << pos << ' '
              << std::setw(10) << tg << ' '

      //<< vol_pre_name << ' ' << vol_postname
              << std::endl;
    return true; // Still don't know what this means!
  };

  n4::sensitive_detector::end_of_event_fn eoe = [&last_interesting_event, &got_electron](auto) {
    if (n4::event_number() == last_interesting_event) {std::cout << std::endl;}
    got_electron = false;
  };

  // pick one:
  auto sd = new n4::sensitive_detector{"Noisy_detector", make_noise, eoe};
  //auto sd = nullptr;

  // Set mandatory initialization classes

  // run_manager takes ownership of geometry
  run_manager -> SetUserInitialization(new n4::geometry{[&phantom, sd]() -> G4VPhysicalVolume* {
    // Pick one (ensure that generator (below) is compatible) ...
    return phantom_in_cylinder(phantom, 200*mm, 40*mm, sd);
    return cylinder_lined_with_hamamatsus(150*mm, 70*mm, 50*mm, sd);
    return phantom.geometry();
    return imas_demonstrator(nullptr);
    return square_array_of_sipms(sd);
    return nain4::place(sipm_hamamatsu_blue(true, sd)).now();
  }});

  { // Physics list
    auto verbosity = 1;
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list -> ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    run_manager  -> SetUserInitialization(physics_list);
  } // run_manager owns physics_list

  // User action initialization
  run_manager->SetUserInitialization(new n4::actions{
      new n4::generator{[&phantom](G4Event* event) {
        // Pick one that matches geometry
        // generate_back_to_back_511_keV_gammas(event, {}, 0);
        phantom.generate_primaries(event);

        // auto optphot = nain4::find_particle("geantino");
        // auto p = G4ThreeVector{0,1,0} * 7 * eV;
        // double time = 0;
        // auto vertex =      new G4PrimaryVertex({0, 65*mm, 3*mm}, time);
        // auto primary = new G4PrimaryParticle(optphot,  p.x(),  p.y(),  p.z());
        // primary -> SetPolarization({1,0,0});
        // vertex->SetPrimary(primary);
        // event -> AddPrimaryVertex(vertex);

      }}});

  // Initialize visualization
  auto vis_manager = make_unique<G4VisExecutive>();
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive{"Quiet"};
  vis_manager -> Initialize();

  // Get the pointer to the User Interface manager
  auto ui_manager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if (!ui) {
    // batch mode
    G4String command = "/control/execute ";
    G4String file_name = argv[1];
    ui_manager -> ApplyCommand(command + file_name);
  } else {
    // interactive mode
    {
      ui_manager -> ApplyCommand("/control/execute init_vis.mac");
      ui_manager -> ApplyCommand("/PhysicsList/RegisterPhysics G4EmStandardPhysics_option4");
      ui_manager -> ApplyCommand("/PhysicsList/RegisterPhysics G4OpticalPhysics");
      ui_manager -> ApplyCommand("/vis/scene/endOfEventAction accumulate 100");
      //ui_manager -> ApplyCommand("/run/beamOn 1");
      //nain4::silence _{G4cout};
      int PHI = 150; int THETA = 160; int THETAF = 175;
      auto view = [&ui_manager](auto theta, auto phi) {
        ui_manager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi "
                                 + std::to_string(theta) + ' ' + std::to_string(phi));
      };
      {
        nain4::silence _{G4cout};
        for (int phi  =PHI  ; phi  <360+PHI   ; phi  +=4) { view(THETA, phi); }
        for (int theta=THETA; theta<360+THETAF; theta+=4) { view(theta, PHI); }
      }
      ui -> SessionStart();
    }
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
}
