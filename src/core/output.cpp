#include "gryphon/core/output.h"

#include "gryphon/utils/io.h"
#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"
#include "gryphon/utils/timer.h"

#define USETHREADS

#ifdef USETHREADS
#include <thread>
#endif

namespace gryphon {
namespace core {

std::string make_filename(std::string initfilename, int id, std::string ext) {
  std::string str = initfilename;
  str += "_";
  str += std::to_string(id);
  str += ".";
  str += ext;
  return str;
}

OutputManager::OutputManager(const Input& input) {
  m_E = utils::LogAxis<double>(input.E_min, input.E_max, input.E_size);
  m_I.resize(input.E_size);
  m_filename = make_filename(input.simname, input.seed, "txt");
}

void OutputManager::dump() const {
  const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out(m_filename);
  out << "# E [GeV] - I [GeV-1 m-2 sec-1 sr-1]\n";
  out << std::scientific;
  for (size_t i = 0; i < m_E.size(); ++i) {
    out << m_E[i] / cgs::GeV << "\t";
    out << m_I[i] / units << "\t";
    out << "\n";
  }
}

void OutputManager::compute(const particle::Particles& particles) {
  utils::Timer timer("running time");
  unsigned int n = m_E.size();
  std::vector<std::thread> threads(n);

#ifdef USETHREADS
  for (unsigned int i = 0; i < n; i++) {
    threads[i] = std::thread(
        [&](int i) {
          double value = 0;
          for (auto& particle : particles) {
            value += particle->get(m_E.at(i));
          }
          m_I.at(i) = value;
        },
        i);
  }

  for (auto& th : threads) {
    th.join();
  }
#else
  for (unsigned int i = 0; i < n; ++i) {
    double value = 0;
    for (auto& particle : particles) {
      value += particle->get(m_E.at(i));
    }
    m_I.at(i) = value;
  }
#endif
}

}  // namespace core
}  // namespace gryphon