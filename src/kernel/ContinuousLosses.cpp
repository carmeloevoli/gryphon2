#include "gryphon/kernel/ContinuousLosses.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace kernel {
namespace {
using Quadrature = std::unique_ptr<gsl_integration_glfixed_table,
                                   decltype(&gsl_integration_glfixed_table_free)>;

Quadrature makeQuadrature(size_t size) {
  Quadrature table(gsl_integration_glfixed_table_alloc(size), gsl_integration_glfixed_table_free);
  if (!table) throw std::runtime_error("Could not allocate steady-source quadrature");
  return table;
}

struct IntegrationPoint { double lambda2; double weight; };

std::vector<IntegrationPoint> integrationPoints(
    double energy, size_t count, const core::Input& input,
    const DiffusionLossesKernel& kernel, const injection::MSPInjection& injection) {
  auto table = makeQuadrature(count);
  const double hi = std::log(40.*input.injEmax()/energy);
  if (!(hi > -32.)) throw std::invalid_argument("Energy/cutoff ratio too large for quadrature");
  std::vector<IntegrationPoint> result;
  result.reserve(count);
  for (size_t j = 0; j < count; ++j) {
    double v, weight;
    gsl_integration_glfixed_point(-32., hi, j, &v, &weight, table.get());
    const double dE = energy*std::exp(v);
    const double es = energy+dE;
    const double lambda2 = kernel.lambda2(energy, es);
    if (!(lambda2 > 0.)) throw std::runtime_error("Invalid diffusion length in steady integral");
    const double normalization = cgs::c_light/(4.*M_PI*kernel.b(energy));
    result.push_back({lambda2, normalization*injection.rate(es, cgs::erg/cgs::second)*
                                  dE*weight/std::pow(M_PI*lambda2, 1.5)});
  }
  return result;
}
}  // namespace

ContinuousLosses::ContinuousLosses(const core::Input& input)
    : m_input(input), m_kernel(input), m_injection(input) {
  input.validate();
  m_energies = utils::LogAxis(input.E_min(), input.E_max(), input.E_size());
  // Exact energies for the source-contribution diagnostics, when in range.
  for (double e : {10., 100., 300., 1000.}) {
    e *= cgs::GeV;
    if (e >= input.E_min() && e <= input.E_max()) {
      const auto near = std::find_if(m_energies.begin(), m_energies.end(), [e](double value) {
        return std::abs(value/e-1.) < 1e-10;
      });
      if (near == m_energies.end()) m_energies.push_back(e);
      else *near = e;
    }
  }
  std::sort(m_energies.begin(), m_energies.end());
  const double lmax = std::sqrt(m_kernel.lambda2(input.E_min(), std::numeric_limits<double>::infinity()));
  const double imageBound = lmax*std::sqrt(-std::log(1e-14))/(2.*input.H());
  if (!std::isfinite(imageBound) || imageBound > 1000.)
    throw std::invalid_argument("Steady image expansion requires too many images; check transport");
  m_images = static_cast<int>(std::ceil(imageBound))+1;
  const double maxDistance = 1.01*std::hypot(input.R_g()+input.R_sun(), (2.*m_images+1.)*input.H());
  m_logDistances = utils::LinAxis(std::log(.1*cgs::pc), std::log(maxDistance), input.mspDistances());
  m_response.resize(m_energies.size(), std::vector<double>(m_logDistances.size()));
  for (size_t e = 0; e < m_energies.size(); ++e) {
    const auto points = integrationPoints(m_energies[e], input.mspQuadrature(), input, m_kernel, m_injection);
    for (size_t d = 0; d < m_logDistances.size(); ++d) {
      const double d2 = std::exp(2.*m_logDistances[d]);
      double value = 0.;
      for (const auto& point : points) value += point.weight*std::exp(-d2/point.lambda2);
      m_response[e][d] = value;
    }
  }
}

double ContinuousLosses::interpolate(size_t energyIndex, double distance) const {
  if (!(distance > 0.)) throw std::invalid_argument("Steady point source coincides with observer");
  const double logd = std::log(distance);
  if (logd < m_logDistances.front() || logd > m_logDistances.back())
    throw std::out_of_range("MSP source/image outside distance table");
  const double fraction = (logd-m_logDistances.front()) /
                          (m_logDistances.back()-m_logDistances.front())*(m_logDistances.size()-1);
  const size_t index = std::min(static_cast<size_t>(fraction), m_logDistances.size()-2);
  const auto& row = m_response.at(energyIndex);
  return row[index]+(fraction-index)*(row[index+1]-row[index]);
}

double ContinuousLosses::response(size_t energyIndex, const utils::Vector3d& pos) const {
  if (std::abs(pos.z) >= m_input.H()) return 0.;
  const double rho2 = pow2(pos.x)+pow2(pos.y);
  const double free = interpolate(energyIndex, pos.getModule());
  double value = free, compensation = 0.;
  for (int n = 1; n <= m_images; ++n) {
    const double sign = n%2 ? -1. : 1.;
    const double dz = 2.*n*m_input.H();
    const double term = sign*(interpolate(energyIndex, std::sqrt(rho2+pow2(sign*pos.z+dz))) +
                              interpolate(energyIndex, std::sqrt(rho2+pow2(sign*pos.z-dz))));
    const double corrected = term-compensation;
    const double next = value+corrected;
    compensation = (next-value)-corrected;
    value = next;
  }
  if (!std::isfinite(value) || value < -1e-8*free)
    throw std::runtime_error("Unstable steady image sum; increase mspdistances");
  return std::max(value, 0.);
}

double ContinuousLosses::directResponse(double energy, const utils::Vector3d& pos, size_t quadrature) const {
  if (std::abs(pos.z) >= m_input.H()) return 0.;
  if (!(pos.getModule() > 0.)) throw std::invalid_argument("Steady point source coincides with observer");
  const auto points = integrationPoints(energy, quadrature, m_input, m_kernel, m_injection);
  double value = 0.;
  for (const auto& point : points) {
    value += point.weight*std::exp(-(pow2(pos.x)+pow2(pos.y))/point.lambda2)*
             utils::halo_function(point.lambda2, m_input.H(), 0., pos.z);
  }
  return value;
}

}  // namespace kernel
}  // namespace gryphon
