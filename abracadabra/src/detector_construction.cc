#include "detector_construction.hh"

#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4Orb.hh>
#include <G4PVPlacement.hh>
#include <G4RunManager.hh>
#include <G4Sphere.hh>
#include <G4String.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Trd.hh>

// Create logical volume from solid and material
G4LogicalVolume* logical(G4Material* material, G4VSolid* solid) {
  return new G4LogicalVolume{solid, material, solid->GetName()};
}

// Utility for concisely creating materials from NIST code
G4Material* material(G4String const& name) { return G4NistManager::Instance()->FindOrBuildMaterial(name); };

// ================================================================================
#include <optional>
using std::nullopt;
using std::optional;
using std::make_optional;

class place {
public:
  place(G4LogicalVolume* child)
  : child(child ? make_optional(child) : nullopt) {}

  place(place const&) = default;

  place& at(G4double x, G4double y, G4double z) { return this->at({x, y, z}); }

  place& at(G4ThreeVector position_) {
    auto xxx = make_optional(position_);
    this->position.swap(xxx);
    return *this;
  }

  place& in(G4LogicalVolume* parent_) {
    auto xxx = make_optional(parent_);
    this->parent.swap(xxx);
    return *this;
  }

  G4PVPlacement* now() {
    // Maybe allow setting these later on
    bool bool_op        = false;
    bool check_overlaps = true;

    return new G4PVPlacement {
      rotation.value_or(nullptr),
      position.value_or(G4ThreeVector{}),
      child   .value(),
      name    .value_or(child.value() -> GetName()),
      parent  .value_or(nullptr),
      bool_op,
      check_overlaps
    };
  }

private:
  optional<G4LogicalVolume*>  child;
  optional<G4LogicalVolume*>  parent;
  optional<G4ThreeVector>     position;
  optional<G4RotationMatrix*> rotation;
  optional<G4String>          name;
};

// ================================================================================

G4VPhysicalVolume* detector_construction::Construct() {

  // ----- Materials --------------------------------------------------------------
  auto air    = material("G4_AIR");
  auto water  = material("G4_WATER");
  auto tissue = material("G4_A-150_TISSUE");
  auto bone   = material("G4_BONE_COMPACT_ICRU");

  // ----- Dimensions -------------------------------------------------------------
  // Size of the detector
  G4double length_xy = 20 * cm;
  G4double length_z  = 30 * cm;

  // Envelope: G4Box requires half-lengths
  G4double e_xy = 0.5 * length_xy;
  G4double e_z  = 0.5 * length_z;

  // World volume needs a margin around everything inside it
  G4double w_xy = 1.2 * e_xy;
  G4double w_z  = 1.2 * e_z;

  // Trapezoid ---------------------------------------------------------------------
  G4double t_dxa = 12 * cm / 2, t_dxb = 12 * cm / 2;
  G4double t_dya = 10 * cm / 2, t_dyb = 16 * cm / 2;
  G4double t_dz  =  6 * cm / 2;

  // Cone --------------------------------------------------------------------------
  G4double c_rmin_a  = 0 * cm,   c_rmax_a = 2 * cm;
  G4double c_rmin_b  = 0 * cm,   c_rmax_b = 4 * cm;
  G4double c_hz      = 3 * cm;
  G4double c_phi_min = 0 * deg,  c_phi_max = 360 * deg;

  // ----- Create the shapes -------------------------------------------------------
  auto world     = logical(air   , new G4Box{"World"   , w_xy, w_xy, w_z});
  auto envelope  = logical(water , new G4Box{"Envelope", e_xy, e_xy, e_z});
  auto trapezoid = logical(bone  , new G4Trd{"BoneTrapezoid", t_dxa, t_dxb, t_dya, t_dyb, t_dz});
  auto cone      = logical(tissue, new G4Cons{"TissueCone",
                                              c_rmin_a, c_rmax_a, c_rmin_b, c_rmax_b,
                                              c_hz, c_phi_min, c_phi_max});

  this->scoring_volume = trapezoid;

  // ----- Place the shapes at specific points in space ----------------------------
  place(trapezoid).in(envelope).at(0, -1*cm, 7*cm).now();
  place(cone     ).in(envelope).at(0,  2*cm,-7*cm).now();
  place(envelope ).in(world   )                   .now();

  return place(world).now();
}
