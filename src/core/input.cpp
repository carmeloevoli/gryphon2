#include "gryphon/core/input.h"

// #include <fstream>
// #include <iostream>
// #include <memory>
// #include <sstream>
// #include <string>

// #include "include/git_revision.h"
// #include "include/utilities.h"

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace core {

Input::Input() {  // validator();
}

Input::Input(const std::string& filename) {
  // read_params_file(filename);
  // validator();
}

// void Input::set_params(const std::string& key, const double& value) {
//   if (key == "emin")
//     _E_min = value * cgs::GeV;
//   else if (key == "emax")
//     _E_max = value * cgs::GeV;
//   else if (key == "esize")
//     _E_size = static_cast<int>(value);
//   else if (key == "hkpc")
//     _H = value * cgs::kpc;
//   else if (key == "rgkpc")
//     _R_g = value * cgs::kpc;
//   else if (key == "sunkpc")
//     _sun = Vector3d(value * cgs::kpc, 0., 0.);
//   else if (key == "va")
//     _v_A = value * cgs::km / cgs::sec;
//   else if (key == "ddelta")
//     _D_delta = value;
//   else if (key == "d0h")
//     _D0_over_H = value * 1e28 * cgs::cm2 / cgs::sec / cgs::kpc;
//   else if (key == "dbreak")
//     _D_E_break = value * cgs::GeV;
//   else if (key == "dhardening")
//     _D_ddelta = value;
//   else if (key == "snralpha")
//     _SNR_alpha = value;
//   else if (key == "snrefficiency")
//     _SNR_efficiency = value;
//   else if (key == "snrcutoff")
//     _SNR_E_c = value * cgs::GeV;
//   else if (key == "meanp0")
//     _pulsar_mean_period = value * cgs::msec;
//   else if (key == "sigmap0")
//     _pulsar_sdev_period = value * cgs::msec;
//   else if (key == "alphalow")
//     _PWN_alpha_low = value;
//   else if (key == "alphahigh")
//     _PWN_alpha_high = value;
//   else if (key == "pwnbreak")
//     _PWN_E_break = value * cgs::GeV;
//   else if (key == "maxtimemyr")
//     _max_time = value * cgs::Myr;
//   else if (key == "rate")
//     _SN_rate = value / (100. * cgs::year);
//   else if (key == "maxradiuskpc")
//     _max_radius = value * cgs::kpc;
//   else if (key == "b0")
//     _b_0 = value * cgs::GeV / cgs::sec;
//   else if (key == "fuv")
//     _f_UV = value;
//   else if (key == "bhalomug" || key == "bmug")
//     _B_halo = value * cgs::microgauss;
//   else
//     throw std::runtime_error("Params " + key + " not found");
// }

// void Input::set_seed(const size_t& seed) { _seed = seed; }

// void eraseSubStr(std::string* main_string, const std::string& sub_string) {
//   // Search for the substring in string
//   //   size_t pos = main_string->find(sub_string);
//   //   if (pos != std::string::npos) {
//   //     main_string->erase(pos, sub_string.length());
//   //   }
// }

// void Input::set_simname(const std::string& params_filename) {
//   //   _simname = params_filename;
//   //   const auto size = _simname.length();
//   //   eraseSubStr(&_simname, "config_files/");
//   //   eraseSubStr(&_simname, ".ini");
//   //   if (_simname.length() != size - 17)
//   //     throw std::runtime_error("Problem with the input filename! Must be :
//   //     config_files/xxx.ini");
// }

// void Input::set_outputdir(const std::string& dir) { _simname = dir + "/" + _simname; }

// void Input::read_params_file(const std::string& params_filename) {
//   //   if (!Utilities::file_exists(params_filename))
//   //     throw std::runtime_error("Input params file does not exist");
//   //   set_simname(params_filename);
//   //   std::ifstream infile(params_filename.c_str());
//   //   std::string line;
//   //   while (std::getline(infile, line)) {
//   //     std::istringstream iss(line);
//   //     std::string key, value;
//   //     if (!(iss >> key >> value)) {
//   //       break;
//   //     }  // error
//   //     key = Utilities::simplify_key(key);
//   //     if (key == "particletype")
//   //       _particleType.set(Utilities::simplify_key(value));
//   //     else if (key == "runmode")
//   //       _runMode.set(Utilities::simplify_key(value));
//   //     else if (key == "lossesmodel")
//   //       _lossesModel.set(Utilities::simplify_key(value));
//   //     else if (key == "spiralmodel")
//   //       _spiralModel.set(Utilities::simplify_key(value));
//   //     else if (key == "snrmodel")
//   //       _snrModel.set(Utilities::simplify_key(value));
//   //     else if (key == "outputdir")
//   //       set_outputdir(value);
//   //     else
//   //       set_params(key, std::stod(value));
//   //   }
// }

// void Input::validator() {
//   //   if (_E_min < cgs::GeV || _E_min > 10.1 * cgs::PeV)
//   //     throw std::runtime_error("E_min outside the range");
//   //   if (_E_max < cgs::GeV || _E_max > 10.1 * cgs::PeV)
//   //     throw std::runtime_error("E_max outside the range");
//   //   if (_E_min >= _E_max) throw std::runtime_error("E_min must be less than E_max");
//   //   if (_E_size < 2) throw std::runtime_error("E_size must be > 1");
//   //   if (_SNR_efficiency > 1 || _SNR_efficiency < 0)
//   //     throw std::runtime_error("SNR efficiency must be < 1");
//   //   if (_SNR_alpha <= 2) throw std::runtime_error("SNR spectrum slope must be > 2.");
//   //   if (_PWN_alpha_low >= 2) throw std::runtime_error("PWN low energy spectrum slope must be
//   //   < 2."); if (_PWN_alpha_high <= 2) throw std::runtime_error("PWN high energy spectrum slope
//   //   must be > 2"); if (_D_E_break < cgs::MeV && _D_E_break > cgs::PeV)
//   //     throw std::runtime_error("Diffusion energy break must be > MeV and < PeV");
//   //   if (_f_UV < 0. || _f_UV > 100.)
//   //     throw std::runtime_error("UV enhancement factor must be > 0 and < 100");
//   //   if (_B_halo < cgs::nanogauss || _B_halo > 10. * cgs::microgauss)
//   //     throw std::runtime_error("B halo must be > nG and < 10 muG");
//   //   // TODO(CE): complete here
// }

// void Input::save_on_file(const std::shared_ptr<RNG> rng) {
//   //   auto filename = Utilities::make_filename(simname, rng->get_seed(), "log");
//   //   std::ofstream outfile(filename.c_str());
//   //   if (outfile.is_open()) {
//   //     std::cout << "Saving input parameters on file " << filename << " ...\n";
//   //     outfile << "git version       : " << Utilities::git_sha1() << "\n";
//   //     outfile << "has local changes : " << std::boolalpha <<
//   Utilities::git_has_local_changes()
//   //             << std::noboolalpha << "\n";
//   //     outfile << "was built on      : " << __DATE__ << " " __TIME__ << "\n";
//   //     outfile << "E_min             : " << _E_min / cgs::GeV << " GeV\n";
//   //     outfile << "E_max             : " << _E_max / cgs::GeV << " GeV\n";
//   //     outfile << "E_size            : " << _E_size << "\n";
//   //     outfile << "halo size         : " << _H / cgs::kpc << " kpc\n";
//   //     outfile << "galaxy radius     : " << _R_g / cgs::kpc << " kpc\n";
//   //     outfile << "sun position      : " << _sun / cgs::kpc << " kpc\n";
//   //     outfile << "v_A               : " << _v_A / (cgs::km / cgs::sec) << " km/s\n";
//   //     outfile << "D delta           : " << _D_delta << "\n";
//   //     outfile << "D break at        : " << _D_E_break / cgs::GeV << " GeV\n";
//   //     outfile << "D ddelta          : " << _D_ddelta << "\n";
//   //     outfile << "D smoothing       : " << _D_smoothness << "\n";
//   //     outfile << "D_0 / H           : " << _D0_over_H / (cgs::cm2 / cgs::sec / cgs::kpc)
//   //             << " cm2/s/kpc\n";
//   //     outfile << "at E_0            : " << _E_0 / cgs::GeV << " GeV\n";
//   //     outfile << "SNR alpha         : " << _SNR_alpha << "\n";
//   //     outfile << "SNR efficiency    : " << _SNR_efficiency << "\n";
//   //     if (_SNR_E_c > 0)
//   //       outfile << "SNR cutoff        : " << _SNR_E_c / cgs::TeV << " TeV\n";
//   //     else
//   //       outfile << "SNR cutoff        : "
//   //               << "none \n";
//   //     outfile << "pulsar_fraction   : " << _pulsar_fraction << "\n";
//   //     outfile << "pulsar mean P0    : " << _pulsar_mean_period / cgs::sec << " s\n";
//   //     outfile << "pulsar sdev P0    : " << _pulsar_sdev_period / cgs::sec << " s\n";
//   //     outfile << "PWN alpha_L       : " << _PWN_alpha_low << "\n";
//   //     outfile << "PWN alpha_H       : " << _PWN_alpha_high << "\n";
//   //     outfile << "PWN E_break       : " << _PWN_E_break / cgs::GeV << " GeV\n";
//   //     outfile << "PWN efficiency    : " << _PWN_efficiency << "\n";
//   //     outfile << "b_0               : " << _b_0 / (cgs::GeV / cgs::sec) << " GeV/s\n";
//   //     outfile << "UV enhancement    : " << _f_UV << "\n";
//   //     outfile << "B halo            : " << _B_halo / cgs::microgauss << "\n";
//   //     outfile << "simulation time   : " << _max_time / cgs::Myr << " Myr\n";
//   //     outfile << "max distance      : " << _max_radius / cgs::kpc << " kpc\n";
//   //     outfile << "SN rate           : " << _SN_rate / (1. / cgs::year) << " yr-1\n";
//   //     outfile << "time step         : " << _time_step / cgs::year << " yr\n";
//   //     outfile << "run mode          : " << _runMode << "\n";
//   //     outfile << "particle type     : " << _particleType << "\n";
//   //     outfile << "losses model      : " << _lossesModel << "\n";
//   //     outfile << "spiral model      : " << _spiralModel << "\n";
//   //     outfile << "SNR model         : " << _snrModel << "\n";
//   //     outfile.close();
//   //     std::cout << "... done.\n";
//   //   } else {
//   //     throw std::runtime_error("parameter file cannot be opened.");
//   //   }
// }

void Input::print() {
  LOGD << "seed : " << _seed;
  LOGD << "E_min : " << _E_min / cgs::GeV << " GeV";
  LOGD << "E_max : " << _E_max / cgs::GeV << " GeV";
  LOGD << "E_size : " << _E_size;
  LOGD << "halo size : " << _H / cgs::kpc << " kpc";
  LOGD << "galaxy radius : " << _Rg / cgs::kpc << " kpc";
  LOGD << "sun position : " << _Rsun / cgs::kpc << " kpc";
  LOGD << "D_0 / H : " << _D0_over_H / (cgs::cm2 / cgs::sec / cgs::kpc) << "cm2/s/kpc";
  LOGD << "E_0 : " << _E_0 / cgs::GeV << " GeV";
  LOGD << "delta : " << _delta;
  LOGD << "ddelta : " << _ddelta;
  LOGD << "s : " << _s;
  LOGD << "E_b : " << _E_b / cgs::GeV << " GeV";
  LOGD << "a : " << _a;
  LOGD << "b : " << _b;
  LOGD << "R_1 : " << _R1 / cgs::kpc << " kpc";
  LOGD << "SN rate : " << _sn_rate / (1. / cgs::year) << " yr-1";
  LOGD << "time step : " << _time_step / cgs::year << " yr";
  LOGD << "max time : " << _max_time / cgs::Myr << " Myr";
}

//   //   std::cout << "galaxy radius     : " << _R_g / cgs::kpc << " kpc\n";
//   //   std::cout << "sun position      : " << _sun / cgs::kpc << " kpc\n";
//   //   std::cout << "SNR alpha         : " << _SNR_alpha << "\n";
//   //   std::cout << "SNR efficiency    : " << _SNR_efficiency << "\n";
//   //   if (_SNR_E_c > 0)
//   //     std::cout << "SNR cutoff        : " << _SNR_E_c / cgs::TeV << " TeV\n";
//   //   else
//   //     std::cout << "SNR cutoff        : "
//   //               << "none \n";
//   //   std::cout << "simulation time   : " << _max_time / cgs::Myr << " Myr\n";
//   //   std::cout << "max distance      : " << _max_radius / cgs::kpc << " kpc\n";
//   //   std::cout << "SN rate           : " << _SN_rate / (1. / cgs::year) << " yr-1\n";
//   //   std::cout << "time step         : " << _time_step / cgs::year << " yr\n";
//   //   std::cout << "run mode          : " << _runMode << "\n";
//   //   std::cout << "particle type     : " << _particleType << "\n";
//   //   std::cout << "losses model      : " << _lossesModel << "\n";
//   //   std::cout << "spiral model      : " << _spiralModel << "\n";
//   //   std::cout << "SNR model         : " << _snrModel << "\n";

}  // namespace core
}  // namespace gryphon