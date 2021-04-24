#include "nain4.hh"

// Solids
#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4Trd.hh>

// Managers
#include <G4NistManager.hh>
#include <G4RunManager.hh>

// Units
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>

// Other G4
#include <G4Material.hh>

#include <catch2/catch.hpp>

// Many of the tests below check physical quantities. Dividing physical
// quantities by their units gives raw numbers which are easily understandable
// by a human reader, which is important test failures are reported. Sometimes
// this gives rise to the apparently superfluous division by the same unit on
// both sides of an equation, in the source code.

#include<numeric>

TEST_CASE("nain material", "[nain][material]") {

  // nain4::material finds the same materials as the verbose G4 style
  SECTION("material NIST") {
    auto material_name = GENERATE("G4_AIR", "G4_WATER", "G4_H", "G4_A-150_TISSUE");
    auto nain_material = nain4::material(material_name);
    auto nist_material = G4NistManager::Instance()->FindOrBuildMaterial(material_name);
    REQUIRE(nain_material == nist_material);
    REQUIRE(nain_material != nullptr);
  }

  // Basic material properties make sense (except for solid water at RTP!)
  SECTION("material properties") {
    SECTION("water") {
      auto water = nain4::material("G4_WATER");
      CHECK(water->GetName()                  == "G4_WATER");
      CHECK(water->GetChemicalFormula()       == "H_2O");
      CHECK(water->GetTemperature() /  kelvin == Approx(293.15));
      CHECK(water->GetPressure() / atmosphere == Approx(1));
      CHECK(water->GetDensity() /     (kg/m3) == Approx(1000));
      CHECK(water->GetState()                 == G4State::kStateSolid); // WTF!?
    }
  }

  // Making and retrieving materials with nain4
  SECTION("material creation from N atoms") {

    // The values used to construct the material
    auto name = "n4test_FR4";
    auto density = 1.85 * g/cm3;
    auto state = kStateSolid;
    auto [nH, nC, nO] = std::make_tuple(12, 18, 3);

    // Make the material using nain4::material_from_elements
    auto fr4 = nain4::material_from_elements_N(name, density, state,
                                               {{"H", nH}, {"C", nC}, {"O", nO}});
    CHECK(fr4 != nullptr);

    // Verify that the material can be retrieved with nani4::material
    auto fr4_found = nain4::material(name);
    CHECK(fr4 == fr4_found);

    // Grab elements and calculate some properties for use in tests lower down
    auto H = nain4::element("H"); auto mH = H->GetAtomicMassAmu();
    auto C = nain4::element("C"); auto mC = C->GetAtomicMassAmu();
    auto O = nain4::element("O"); auto mO = O->GetAtomicMassAmu();
    auto total_mass = nH*mH + nC*mC + nO*mO;

    // Elements used correctly?
    CHECK(fr4 -> GetElement(0) == H);
    CHECK(fr4 -> GetElement(1) == C);
    CHECK(fr4 -> GetElement(2) == O);

    // Correct number of each element?
    auto atoms = fr4 -> GetAtomsVector();
    CHECK(atoms[0] == nH);
    CHECK(atoms[1] == nC);
    CHECK(atoms[2] == nO);

    // Basic properties set corretly?
    CHECK(fr4 -> GetNumberOfElements() == 3);
    CHECK(fr4 -> GetDensity() == density);
    CHECK(fr4 -> GetState() == state);

    // Fractional composition correct?
    auto fracs = fr4 -> GetFractionVector();
    CHECK(fracs[0] == Approx(nH*mH / total_mass));
    CHECK(fracs[1] == Approx(nC*mC / total_mass));
    CHECK(fracs[2] == Approx(nO*mO / total_mass));

    // Does fractional composition sum to 1?
    CHECK(std::accumulate(fracs, fracs + fr4->GetNumberOfElements(), 0.0) == Approx(1));
  }

  // Making and retrieving materials with nain4
  SECTION("material creation from mass fractions") {

    // The values used to construct the material
    auto name = "n4test_LYSO";
    auto density = 7.1 * g/cm3;
    auto state = kStateSolid;
    auto [fLu, fY, fSi, fO] = std::make_tuple(0.714, 0.040, 0.064, 0.182);

    // Make the material using nain4::material_from_elements
    auto lyso = nain4::material_from_elements_F(name, density, state,
                                                {{"Lu", fLu}, {"Y", fY}, {"Si", fSi}, {"O", fO}});
    CHECK(lyso != nullptr);

    // Verify that the material can be retrieved with nani4::material
    auto fr4_found = nain4::material(name);
    CHECK(lyso == fr4_found);

    // Grab elements and calculate some properties for use in tests lower down
    auto Lu = nain4::element("Lu");
    auto Y  = nain4::element("Y" );
    auto Si = nain4::element("Si");
    auto O  = nain4::element("O" );
    //auto total_mass = nH*mH + nC*mC + nO*mO;

    // Elements used correctly?
    CHECK(lyso -> GetElement(0) == Lu);
    CHECK(lyso -> GetElement(1) == Y );
    CHECK(lyso -> GetElement(2) == Si);
    CHECK(lyso -> GetElement(3) == O );

    // Atom counts produce nonsense when material built with mass fractions
    auto atoms = lyso -> GetAtomsVector();
    CHECK(atoms[0] == 1);
    CHECK(atoms[1] == 0);
    CHECK(atoms[2] == 0);
    CHECK(atoms[3] == 2);

    // Basic properties set corretly?
    CHECK(lyso -> GetNumberOfElements() == 4);
    CHECK(lyso -> GetDensity() == density);
    CHECK(lyso -> GetState() == state);

    // Fractional composition correct?
    auto fracs = lyso -> GetFractionVector();
    CHECK(fracs[0] == fLu);
    CHECK(fracs[1] == fY );
    CHECK(fracs[2] == fSi);
    CHECK(fracs[3] == fO );

    // Does fractional composition sum to 1?
    CHECK(std::accumulate(fracs, fracs + lyso->GetNumberOfElements(), 0.0) == Approx(1));

  }
}

TEST_CASE("nain volume", "[nain][volume]") {
  // nain4::volume produces objects with sensible sizes, masses, etc.
  auto water = nain4::material("G4_WATER");
  auto lx = 1 * m;
  auto ly = 2 * m;
  auto lz = 3 * m;
  auto box = nain4::volume<G4Box>("test_box", water, lx, ly, lz);
  auto density = water->GetDensity();
  CHECK(box->TotalVolumeEntities() == 1);
  CHECK(box->GetMass() / kg        == Approx(8 * lx * ly * lz * density / kg));
  CHECK(box->GetMaterial()         == water);
  CHECK(box->GetName()             == "test_box");

  auto solid = box->GetSolid();
  CHECK(solid->GetCubicVolume() / m3 == Approx(8 *  lx    * ly    * lz     / m3));
  CHECK(solid->GetSurfaceArea() / m2 == Approx(8 * (lx*ly + ly*lz + lz*lx) / m2));
  CHECK(solid->GetName()             == "test_box");
}

TEST_CASE("nain place", "[nain][place]") {
  // nain4::place is a good replacement for G4PVPlacement
  auto air = nain4::material("G4_AIR");
  auto outer = nain4::volume<G4Box>("outer", air, 1*m, 2*m, 3*m);

  // Default values are sensible
  SECTION("defaults") {
    auto world = nain4::place(outer).now();

    auto trans = world->GetObjectTranslation();
    CHECK(trans == G4ThreeVector{});
    CHECK(world -> GetName()          == "outer");
    CHECK(world -> GetCopyNo()        == 0);
    CHECK(world -> GetLogicalVolume() == outer);
    CHECK(world -> GetMotherLogical() == nullptr);
  }

  // Multiple optional values can be set at once.
  SECTION("multiple options") {
    G4ThreeVector translation = {1,2,3};
    auto world = nain4::place(outer)
      .at(translation) // 1-arg version of at()
      .name("not outer")
      .copy_no(382)
      .now();

    CHECK(world -> GetObjectTranslation() == translation);
    CHECK(world -> GetName()   == "not outer");
    CHECK(world -> GetCopyNo() == 382);
  }

  // The at() option accepts vector components (as well as a whole vector)
  SECTION("at 3-args") {
    auto world = nain4::place(outer).at(4,5,6).now(); // 3-arg version of at()
    CHECK(world->GetObjectTranslation() == G4ThreeVector{4,5,6});
  }

  // The in() option creates correct mother/daughter relationship
  SECTION("in") {
    auto water = nain4::material("G4_WATER");
    auto inner = nain4::volume<G4Box>("inner", water, 0.3*m, 0.2*m, 0.1*m);

    auto inner_placed = nain4::place(inner)
      .in(outer)
      .at(0.1*m, 0.2*m, 0.3*m)
      .now();

    auto outer_placed = nain4::place(outer).now();

    CHECK(inner_placed -> GetMotherLogical() == outer);
    CHECK(outer_placed -> GetLogicalVolume() == outer);
    CHECK(outer -> GetNoDaughters() == 1);
    CHECK(outer -> GetDaughter(0) -> GetLogicalVolume() == inner);

    // Quick visual check that geometry_iterator works TODO expand
    SECTION("geometry iterator") {
      std::cout << std::endl;
      for (const auto v: outer_placed) {
        std::cout << std::setw(15) << v->GetName() << ": ";
        auto l = v->GetLogicalVolume();
        std::cout
          << std::setw(12) << l->GetMaterial()->GetName()
          << std::setw(12) << G4BestUnit(l->GetMass(), "Mass")
          << std::setw(12) << G4BestUnit(l->GetSolid()->GetCubicVolume(), "Volume")
          << std::endl;
      }
      std::cout << std::endl;
    }
  }
}

TEST_CASE("nain scale_by", "[nain][scale_by]") {
  CHECK(nain4::scale_by(eV, {1, 2.3, 4.5}) == std::vector<G4double>{1*eV, 2.3*eV, 4.5*eV});
  CHECK(nain4::scale_by(cm, {6, 7})        == std::vector<G4double>{6*cm, 7*cm});
}

TEST_CASE("nain vis_attributes", "[nain][vis_attributes]") {
  // Utility for more convenient configuration of G4VisAttributes
  // TODO could do with more extensive testing
  using nain4::vis_attributes;
  auto convenient = vis_attributes{}
    .visible(true)
    .colour({1,0,0})
    .start_time(1.23)
    .end_time(4.56)
    .force_line_segments_per_circle(20)
    .force_solid(true)
    .force_wireframe(false);
  auto pita = G4VisAttributes{};
  pita.SetVisibility(true);
  pita.SetColour({1,0,0});
  pita.SetStartTime(1.23);
  pita.SetEndTime(4.56);
  pita.SetForceLineSegmentsPerCircle(20);
  pita.SetForceSolid(true);
  pita.SetForceWireframe(false);
  CHECK(convenient == pita);

  // The meaning of the different constructors

  // Default constructor sets colour: white, visibility: true
  CHECK(vis_attributes{} == vis_attributes{}.colour({1,1,1}).visible(true));
  // Can set colour via constructor
  CHECK(vis_attributes         {{0,1,0}} ==
        vis_attributes{}.colour({0,1,0}));
  // Can set visibility via constructor
  CHECK(vis_attributes          {true} ==
        vis_attributes{}.visible(true));

  CHECK(vis_attributes          {false} ==
        vis_attributes{}.visible(false));
  // Can set both visibility and colour via constructor
  CHECK(vis_attributes          {false ,       {1,1,0}} ==
        vis_attributes{}.visible(false).colour({1,1,0}));
}

TEST_CASE("nain find", "[nain][find]") {
  // Utilities for retrieving from stores
  SECTION("find_logical") {
    auto air       = nain4::material("G4_AIR");
    auto long_name = "made just for find_logical test";
    auto find_me   = nain4::volume<G4Box>(long_name, air, 1 * cm, 1 * cm, 1 * cm);
    auto found     = nain4::find_logical(long_name);
    CHECK(found == find_me);
    auto should_not_exist = nain4::find_logical("Hopefully this name hasn't been used anywhere", false);
    CHECK(should_not_exist == nullptr);

    SECTION("find_physical") {
      auto placed = nain4::place(find_me).now();
      auto found_placed = nain4::find_physical(long_name);
      CHECK(found_placed == placed);
      CHECK(found_placed != nullptr);
    }
  }

  SECTION("find_particle") {
    auto name = "gamma";
    auto pita = G4ParticleTable::GetParticleTable()->FindParticle(name);
    auto convenient = nain4::find_particle(name);
    CHECK(convenient == pita);
  }
}

TEST_CASE("nain clear_geometry", "[nain][clear_geometry]") {
  auto name = "vanish";
  auto air = nain4::material("G4_AIR");
  auto logical = nain4::volume<G4Box>(name, air, 1*cm, 1*cm, 1*cm);
  auto solid = logical -> GetSolid();
  auto physical = nain4::place(logical).now();
  auto verbose = false;
  {
    auto found_solid    = nain4::find_solid   (name, verbose);
    auto found_logical  = nain4::find_logical (name, verbose);
    auto found_physical = nain4::find_physical(name, verbose);
    CHECK(found_solid    != nullptr);
    CHECK(found_logical  != nullptr);
    CHECK(found_physical != nullptr);
    CHECK(found_solid    == solid);
    CHECK(found_logical  == logical);
    CHECK(found_physical == physical);
  }
  // Clear geometry and verify they are all gone
  nain4::clear_geometry();
  {
    auto found_solid    = nain4::find_solid   (name, verbose);
    auto found_logical  = nain4::find_logical (name, verbose);
    auto found_physical = nain4::find_physical(name, verbose);
    CHECK(found_solid    == nullptr);
    CHECK(found_logical  == nullptr);
    CHECK(found_physical == nullptr);
  }

}

TEST_CASE("nain geometry iterator", "[nain][geometry][iterator]") {
  auto air = nain4::material("G4_AIR");

  auto l   = nain4::volume<G4Box>("l",   air, 100*m, 100*m, 100*m);
  auto l1  = nain4::volume<G4Box>("l1",  air,  40*m,  40*m,  40*m);
  auto l2  = nain4::volume<G4Box>("l2",  air,  10*m,  10*m,  10*m);
  auto l11 = nain4::volume<G4Box>("l11", air,  10*m,  10*m,  10*m);
  auto l21 = nain4::volume<G4Box>("l21", air,  10*m,  10*m,  10*m);
  auto l22 = nain4::volume<G4Box>("l22", air,  10*m,  10*m,  10*m);

  auto p   = nain4::place(l  )                           .now();
  auto p1  = nain4::place(l1 ).in(l) .at(-30*m, 0*m, 0*m).now();
  auto p2  = nain4::place(l2 ).in(l) .at( 30*m, 0*m, 0*m).now();
  auto p11 = nain4::place(l11).in(l1)                    .now();
  auto p21 = nain4::place(l21).in(l2).at(-20*m, 0*m, 0*m).now();
  auto p22 = nain4::place(l22).in(l2).at( 20*m, 0*m, 0*m).now();

  std::vector<G4VPhysicalVolume*> found{begin(p), end(p)};
  std::vector<G4VPhysicalVolume*> expected{p, p1, p2, p11, p21, p22};
  CHECK(found == expected);

}

#include <G4UImanager.hh>

TEST_CASE("nain basic messenger", "[nain][messenger]") {
  GIVEN ("the simplest, unitless messenger") {
    class Dummy{
    public:
      int         foo;
      double      bar;
      std::string baz;
      bool        qux;
      std::unique_ptr<nain4::Messenger> msg;

      Dummy() : foo(0), bar(0), baz(""), qux(false), msg(nullptr) {
        msg = std::unique_ptr<nain4::Messenger>{new nain4::Messenger{this, "/dummy/", "msg doc"}};
        msg->add("foo", foo, "foo doc");
        msg->add("bar", bar, "bar doc");
        msg->add("baz", baz, "baz doc");
        msg->add("qux", qux, "qux doc");
      }

    };
  auto dummy = Dummy{};

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/dummy/foo 666");
  UI->ApplyCommand("/dummy/bar 3.1416");
  UI->ApplyCommand("/dummy/baz spaghetti");
  UI->ApplyCommand("/dummy/qux true");

  CHECK(dummy.foo ==           666 );
  CHECK(dummy.bar == Approx(3.1416));
  CHECK(dummy.baz ==    "spaghetti");
  CHECK(dummy.qux                  );
  }



  GIVEN ("a messenger with units") {
    class Dummy{
    public:
      double foo;
      double bar;
      double baz;
      std::unique_ptr<nain4::Messenger> msg;

      Dummy() : foo(0), bar(0), baz(0), msg(nullptr) {
        msg = std::unique_ptr<nain4::Messenger>{new nain4::Messenger{this, "/dummy/", "msg doc"}};
        msg->add("foo", foo, "foo doc").unit("Length");
        msg->add("bar", bar, "bar doc").unit("Energy");
        msg->add("baz", baz, "baz doc").unit("Time"  );
      }

    };
  auto dummy = Dummy{};

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/dummy/foo 1.5708  cm");
  UI->ApplyCommand("/dummy/bar 3.1416 TeV");
  UI->ApplyCommand("/dummy/baz 6.2832  ps");

  CHECK(dummy.foo == Approx(1.5708 *  cm));
  CHECK(dummy.bar == Approx(3.1416 * TeV));
  CHECK(dummy.baz == Approx(6.2832 *  ps));
  }


  GIVEN ("a messenger with auto units") {
    class Dummy{
    public:
      double foo;
      double bar;
      double baz;
      double qux;
      double quux;
      double corge;
      double grault;
      double garply;
      double waldo;

      std::unique_ptr<nain4::Messenger> msg;

      Dummy() : msg(nullptr) {
        msg = std::unique_ptr<nain4::Messenger>{new nain4::Messenger{this, "/dummy/", "msg doc"}};
        msg->add("radius"                  ,    foo);
        msg->add("bar_diameter"            ,    bar);
        msg->add("the_length_of_baz"       ,    baz);
        msg->add("qux_height"              ,    qux);
        msg->add("the_sides_of_the_objects",   quux);
        msg->add("some_wavelengths"        ,  corge);
        msg->add("thickness_of_a_volume"   , grault);

        msg->add("garply_energy"           , garply);

        msg->add("time_waldo"              ,  waldo);

      }

    };
  auto dummy = Dummy{};

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/dummy/radius                   3.1416 cm");
  UI->ApplyCommand("/dummy/bar_diameter             3.1416 cm");
  UI->ApplyCommand("/dummy/the_length_of_baz        3.1416 cm");
  UI->ApplyCommand("/dummy/qux_height               3.1416 cm");
  UI->ApplyCommand("/dummy/the_sides_of_the_objects 3.1416 cm");
  UI->ApplyCommand("/dummy/some_wavelengths         3.1416 cm");
  UI->ApplyCommand("/dummy/thickness_of_a_volume    3.1416 cm");

  UI->ApplyCommand("/dummy/garply_energy 6.2832 GeV");

  UI->ApplyCommand("/dummy/time_waldo 1.5708 ps");

  CHECK(dummy.foo    == Approx(3.1416 * cm));
  CHECK(dummy.bar    == Approx(3.1416 * cm));
  CHECK(dummy.baz    == Approx(3.1416 * cm));
  CHECK(dummy.qux    == Approx(3.1416 * cm));
  CHECK(dummy.quux   == Approx(3.1416 * cm));
  CHECK(dummy.corge  == Approx(3.1416 * cm));
  CHECK(dummy.grault == Approx(3.1416 * cm));
  CHECK(dummy.garply == Approx(6.2832 * GeV));
  CHECK(dummy.waldo  == Approx(1.5708 *  ps));
  }


  GIVEN ("a messenger with restricted range without units and correct values ") {
    class Dummy{
    public:
      double foo;
      double bar;
      double baz;
      double qux;
      double quux;
      double corge;
      double grault;
      double garply;
      double waldo;
      double fred;
      double plugh;

      std::unique_ptr<nain4::Messenger> msg;

      Dummy() : msg(nullptr) {
        msg = std::unique_ptr<nain4::Messenger>{new nain4::Messenger{this, "/dummy/", "msg doc"}};
        msg->add(   "foo",    foo,    "foo doc").gt      (1);
        msg->add(   "bar",    bar,    "bar doc").gteq    (2);
        msg->add(   "baz",    baz,    "baz doc").lt      (3);
        msg->add(   "qux",    qux,    "qux doc").lteq    (4);
        msg->add(  "quux",   quux,   "quux doc").range   (5,  7);
        msg->add( "corge",  corge,  "corge doc").exclude (8, 10);
        msg->add("grault", grault, "grault doc").positive();
        msg->add("garply", garply, "garply doc").strictly_positive();
        msg->add( "waldo",  waldo,  "waldo doc").negative();
        msg->add(  "fred",   fred,   "fred doc").strictly_negative();
        // undefined reference to `Command& Command::one_of<double>(std::initializer_list<double>)'
        // WTF
        // msg->add( "plugh",  plugh,  "plugh doc").unit("Time"  ).one_of  ({1*s, 5*s, 9*s});
      }

    };
  auto dummy = Dummy{};

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/dummy/foo     2");
  UI->ApplyCommand("/dummy/bar     2");
  UI->ApplyCommand("/dummy/baz     2");
  UI->ApplyCommand("/dummy/qux     4");
  UI->ApplyCommand("/dummy/quux    6");
  UI->ApplyCommand("/dummy/corge   5");
  UI->ApplyCommand("/dummy/grault  0");
  UI->ApplyCommand("/dummy/garply  1");
  UI->ApplyCommand("/dummy/waldo   0");
  UI->ApplyCommand("/dummy/fred   -1");
  // UI->ApplyCommand("/dummy/plugh   9   s");

  CHECK(dummy.foo    ==  2);
  CHECK(dummy.bar    ==  2);
  CHECK(dummy.baz    ==  2);
  CHECK(dummy.qux    ==  4);
  CHECK(dummy.quux   ==  6);
  CHECK(dummy.corge  ==  5);
  CHECK(dummy.grault ==  0);
  CHECK(dummy.garply ==  1);
  CHECK(dummy.waldo  ==  0);
  CHECK(dummy.fred   == -1);
  // CHECK(dummy.plugh  == Approx( 9 *   s));
  }

  GIVEN ("a messenger with restricted range with units and correct values") {
    class Dummy{
    public:
      double foo;
      double bar;
      double baz;
      double qux;
      double quux;
      double corge;
      double grault;
      double garply;
      double waldo;
      double fred;
      int    plugh;

      std::unique_ptr<nain4::Messenger> msg;

      Dummy() : msg(nullptr) {
        msg = std::unique_ptr<nain4::Messenger>{new nain4::Messenger{this, "/dummy/", "msg doc"}};
        msg->add(   "foo",    foo,    "foo doc").unit("Length").gt      (1 * cm);
        msg->add(   "bar",    bar,    "bar doc").unit("Length").gteq    (2 * cm);
        msg->add(   "baz",    baz,    "baz doc").unit("Energy").lt      (3 * keV);
        msg->add(   "qux",    qux,    "qux doc").unit("Energy").lteq    (4 * keV);
        msg->add(  "quux",   quux,   "quux doc").unit("Time"  ).range   (1 * ns, 1 * us);
        msg->add( "corge",  corge,  "corge doc").unit("Time"  ).exclude (1 * ns, 1 * us);
        msg->add("grault", grault, "grault doc").unit("Length").positive();
        msg->add("garply", garply, "garply doc").unit("Length").strictly_positive();
        msg->add( "waldo",  waldo,  "waldo doc").unit("Energy").negative();
        msg->add(  "fred",   fred,   "fred doc").unit("Energy").strictly_negative();
        // undefined reference to `Command& Command::one_of<double>(std::initializer_list<double>)'
        // WTF
        // msg->add( "plugh",  plugh,  "plugh doc").unit("Time"  ).one_of  ({1*s, 5*s, 9*s});
      }

    };
  auto dummy = Dummy{};

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/dummy/foo     1   m");
  UI->ApplyCommand("/dummy/bar     2  cm");
  UI->ApplyCommand("/dummy/baz     1  eV");
  UI->ApplyCommand("/dummy/qux     4 keV");
  UI->ApplyCommand("/dummy/quux    5  ns");
  UI->ApplyCommand("/dummy/corge   1  ms");
  UI->ApplyCommand("/dummy/grault  0  nm");
  UI->ApplyCommand("/dummy/garply  1  nm");
  UI->ApplyCommand("/dummy/waldo   0  eV");
  UI->ApplyCommand("/dummy/fred   -1 keV");
  // UI->ApplyCommand("/dummy/plugh   9   s");

  CHECK(dummy.foo    == Approx( 1 *   m));
  CHECK(dummy.bar    == Approx( 2 *  cm));
  CHECK(dummy.baz    == Approx( 1 *  eV));
  CHECK(dummy.qux    == Approx( 4 * keV));
  CHECK(dummy.quux   == Approx( 5 *  ns));
  CHECK(dummy.corge  == Approx( 1 *  ms));
  CHECK(dummy.grault == Approx( 0 *  nm));
  CHECK(dummy.garply == Approx( 1 *  nm));
  CHECK(dummy.waldo  == Approx( 0 *  eV));
  CHECK(dummy.fred   == Approx(-1 * keV));
  // CHECK(dummy.plugh  == Approx( 9 *   s));
  }



}

/*
TEST_CASE("nain messenger", "[nain][messenger]") {
  GIVEN ("a dummy class that uses a messenger") {
    class Dummy{
    public:
      double fFoo;
      Messenger fmsg;

      Dummy(){
        fmsg = Messenger{this, "/path", "msgdoc"};
        fmsg.add("a_command", fFoo, "cmddoc");
      }
    };

    THEN("this class should throw an exception when initialized"){
      CHECK_THROWS(Dummy());
    }
  }
}

*/
