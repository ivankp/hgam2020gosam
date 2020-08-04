#include <chrono>
#include <cmath>

#include "Higgs2diphoton.hh"

Higgs2diphoton::Higgs2diphoton(seed_type seed)
: gen(seed ?: std::chrono::system_clock::now().time_since_epoch().count()),
  phi_dist(0.,2*M_PI), cts_dist(-1.,1.)
{ }

Higgs2diphoton::photons_type
Higgs2diphoton::operator()(const vec_t& Higgs, bool new_kin) {
  if (new_kin) {
    const double phi = phi_dist(gen);
    const double cts = cts_dist(gen);

    const double sts = std::sin(std::acos(cts));
    const double cos_phi = std::cos(phi);
    const double sin_phi = std::sin(phi);

    photon = { cos_phi*sts, sin_phi*sts, cts };
  }

  const double E = Higgs.m()/2;
  const auto boost = Higgs.boost_vector();

  auto new_photon = photon;
  new_photon.rotate_u_z(boost.normalized()) *= E;

  photons_type diphoton {{ {photon,E}, {-photon,E} }};
  std::get<0>(diphoton) >> boost;
  std::get<1>(diphoton) >> boost;

  return diphoton;
}
