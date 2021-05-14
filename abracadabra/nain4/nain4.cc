#include "nain4.hh"

#include <G4EmStandardPhysics_option4.hh>
#include <G4OpticalPhysics.hh>
#include <FTFP_BERT.hh>
#include <G4PVPlacement.hh>
#include <G4String.hh>

#include <algorithm>
#include <iterator>

namespace nain4 {

G4PVPlacement* place::now() {
  // ----- Name --------------------------------------------------
  // + By default, the name is copied from the child volume.
  // + If a copy_number is specified, it is appended to the name.
  // + All of this is overriden if a name is provided explicitly.
  G4String the_name;
  if (this->label) {
    the_name = this->label.value();
  } else {
    the_name = this->child.value()->GetName();
    // TODO: G4 already appends the copy number to the name?
    if (this->copy_number) {
      auto suffix = "-" + std::to_string(copy_number.value());
      the_name += suffix;
    }
  }
  // TODO: Think about these later
  bool WTF_is_pMany   = false;
  bool check_overlaps = false;

  return new G4PVPlacement{transformation,
                           child.value(),
                           the_name,
                           parent.value_or(nullptr),
                           WTF_is_pMany,
                           copy_number.value_or(0),
                           check_overlaps};
}

std::vector<G4double> scale_by(G4double factor, std::initializer_list<G4double> const& data) {
  std::vector<G4double> out;
  out.reserve(data.size());
  std::transform(begin(data), end(data), back_inserter(out), [factor](auto d){ return d*factor; });
  return out;
}

// --------------------------------------------------------------------------------
// definition of material_properties
material_properties& material_properties::add(G4String const& key, vec const& energies, vec const& values) {
  table -> AddProperty(key, energies, values); // es-vs size equality assertion done in AddProperty
  return *this;
}

material_properties& material_properties::add(G4String const& key, vec const& energies, G4double   value ) {
  return add(key, energies, vec(energies.size(), value));
}

material_properties& material_properties::add(G4String const& key, G4double value) {
  table->AddConstProperty(key, value);
  return *this;
}

// --------------------------------------------------------------------------------
// definition of sensitive_detector
sensitive_detector::sensitive_detector(G4String name, process_hits_fn process_hits, end_of_event_fn end_of_event)
: G4VSensitiveDetector{name}
, process_hits{process_hits}
, end_of_event{end_of_event} {
  fully_activate_sensitive_detector(this);
};

G4bool sensitive_detector::ProcessHits(G4Step* step, G4TouchableHistory*) {
  return process_hits ? process_hits(step) : true;  // TODO what does true mean?
}

void sensitive_detector::EndOfEvent (G4HCofThisEvent* hc) {
  if (end_of_event) { end_of_event(hc); }
}

// --------------------------------------------------------------------------------
// stream redirection utilities

// redirect to arbitrary stream or buffer
redirect::redirect(std::ios& stream, std::streambuf* new_buffer)
: original_buffer(stream.rdbuf())
, stream(stream) {
  stream.rdbuf(new_buffer);
}

redirect::redirect(std::ios& stream, std::ios& new_stream) : redirect{stream, new_stream.rdbuf()} {}

redirect::~redirect() { stream.rdbuf(original_buffer); }

// redirect to /dev/null
silence::silence(std::ios& stream)
  : original_buffer{stream.rdbuf()}
  , stream{stream}
  , dev_null{"/dev/null"} {
  stream.rdbuf(dev_null.rdbuf());
}
silence::~silence() { stream.rdbuf(original_buffer); }


// Set optical physics lists
void use_our_optical_physics(G4RunManager* run_manager, G4int verbosity) {
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list -> ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    run_manager  -> SetUserInitialization(physics_list);
} // run_manager owns physics_list

// --------------------------------------------------------------------------------
// definition of Messenger
Command Messenger::add( const G4String & name
                      , const G4AnyType& var
                      , const G4String & doc // = ""
                      ){
  auto cmd = Command{messenger.DeclareProperty(name, var, doc), name};
  commands.emplace(name, cmd);
  return cmd;
}


Command::Command(G4GenericMessenger::Command& cmd, const G4String cname) :
command(cmd), name(cname) {
  command.SetParameterName(cname, false);
  auto_unit();
}

template<typename T>
Command& Command::one_of(const std::initializer_list<T> values){
  G4String condition = "";
  for (auto v : values){
    condition += compare("==", v) + "||";
  }
  // TODO: Ugly. And unsafe for empty vals
  condition.pop_back();
  condition.pop_back();

  command.SetRange(condition);
  return *this;
}


#define length_keywords {"distance", "radius", "diameter", "height", "side", "length", "wavelength", "thickness"}
#define energy_keywords {"energy"}
#define   time_keywords {"time"}

void Command::auto_unit(){
  for (auto kw : length_keywords){ if (name.contains(kw)) { unit("Length"); return; } }
  for (auto kw : energy_keywords){ if (name.contains(kw)) { unit("Energy"); return; } }
  for (auto kw :   time_keywords){ if (name.contains(kw)) { unit("Time"  ); return; } }
}

#undef length_keywords
#undef energy_keywords
#undef   time_keywords


} // namespace nain4

geometry_iterator begin(G4VPhysicalVolume& vol) { return geometry_iterator{&vol}; }
geometry_iterator   end(G4VPhysicalVolume&    ) { return geometry_iterator{    }; }
geometry_iterator begin(G4VPhysicalVolume* vol) { return begin(*vol); }
geometry_iterator   end(G4VPhysicalVolume* vol) { return   end(*vol); }

geometry_iterator begin(G4LogicalVolume& vol) { return geometry_iterator{&vol}; }
geometry_iterator   end(G4LogicalVolume&    ) { return geometry_iterator{    }; }
geometry_iterator begin(G4LogicalVolume* vol) { return begin(*vol); }
geometry_iterator   end(G4LogicalVolume* vol) { return   end(*vol); }



///////////////////////////////////////////////
