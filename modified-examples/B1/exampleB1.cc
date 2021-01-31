#include "detector_construction.hh"
#include "action_initialization.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "QBBC.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include <G4RunManager.hh>
#include <G4VisManager.hh>

#include <memory>

using std::unique_ptr;

int main(int argc, char **argv) {
  // Detect interactive mode (if no arguments) and define UI session
  auto ui = unique_ptr<G4UIExecutive>{};
  if ( argc == 1 ) {
    ui.reset(new G4UIExecutive(argc, argv));
  }

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);

  // Construct the default run manager
  auto runManager = unique_ptr<G4RunManager>
    {G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default)};

  // Set mandatory initialization classes

  // Detector construction
  runManager -> SetUserInitialization(new detector_construction());

  // Physics list
  G4VModularPhysicsList* physicsList = new QBBC;
  physicsList -> SetVerboseLevel(1);
  runManager  -> SetUserInitialization(physicsList);

  // User action initialization
  runManager -> SetUserInitialization(new action_initialization());

  // Initialize visualization
  auto visManager = unique_ptr<G4VisManager>{new G4VisExecutive};
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager -> Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if (!ui) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager -> ApplyCommand(command+fileName);
  } else {
    // interactive mode
    UImanager -> ApplyCommand("/control/execute init_vis.mac");
    ui -> SessionStart();
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
}
