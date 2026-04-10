#include <unistd.h>

#include <cstdio>
#include <fstream>
#include <string>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {
namespace {

std::string makeTempFilePath() {
  char path[] = "/tmp/gryphon_input_XXXXXX";
  const int fd = mkstemp(path);
  if (fd < 0) throw std::runtime_error("failed to create temporary input file");
  close(fd);
  return std::string(path);
}

}  // namespace

TEST(InputValidation, DefaultConfigurationIsValid) {
  core::Input in;
  EXPECT_NO_THROW(in.validate());
}

TEST(InputValidation, ReportsMultipleErrorsInOneMessage) {
  core::Input in;
  in.set_simname("");
  in.set_simEmin(10. * cgs::TeV);
  in.set_simEmax(1. * cgs::TeV);
  in.set_simEsize(1);
  in.set_rate(2.0 / cgs::year);

  try {
    in.validate();
    FAIL() << "Expected std::invalid_argument";
  } catch (const std::invalid_argument& ex) {
    const std::string message = ex.what();
    EXPECT_NE(message.find("simname cannot be empty"), std::string::npos);
    EXPECT_NE(message.find("E_max must be greater than E_min"), std::string::npos);
    EXPECT_NE(message.find("E_size must be >= 2"), std::string::npos);
    EXPECT_NE(message.find("time_step * sn_rate must be <= 1 for galaxy generation"),
              std::string::npos);
  }
}

TEST(InputValidation, RejectsSunOutsideGalaxy) {
  core::Input in;
  in.set_sunRadius(25. * cgs::kpc);
  EXPECT_THROW(in.validate(), std::invalid_argument);
}

TEST(InputValidation, RejectsSimulationEnergyFloorAtOrBelowOneGeV) {
  core::Input in;
  in.set_simEmin(1.0 * cgs::GeV);
  EXPECT_THROW(in.validate(), std::invalid_argument);
}

TEST(InputValidation, RejectsInvalidInjectionCutoff) {
  core::Input in;
  in.set_injectionModel(InjectionModel::SinglePowerLaw);
  in.set_refEnergy(1.0 * cgs::TeV);
  in.set_injEmax(0.5 * cgs::GeV);
  EXPECT_THROW(in.validate(), std::invalid_argument);
}

TEST(InputValidation, AcceptsInjectionCutoffBelowReferenceEnergyWhenAboveFixedThreshold) {
  core::Input in;
  in.set_injectionModel(InjectionModel::SinglePowerLaw);
  in.set_refEnergy(1.0 * cgs::TeV);
  in.set_injEmax(10.0 * cgs::GeV);
  EXPECT_NO_THROW(in.validate());
}

TEST(InputValidation, RejectsInvalidRandomEmaxBounds) {
  core::Input in;
  in.set_injectionModel(InjectionModel::RandomEmax);
  in.set_injEmax(200.0 * cgs::TeV);
  in.set_injEmaxMin(300.0 * cgs::TeV);
  in.set_injEmaxMax(100.0 * cgs::TeV);
  EXPECT_THROW(in.validate(), std::invalid_argument);
}

TEST(InputParsing, ReadsTextConfigurationFile) {
  const auto path = makeTempFilePath();
  {
    std::ofstream file(path);
    ASSERT_TRUE(file.is_open());
    file << "# Sample config\n";
    file << "simname = parsed_run\n";
    file << "seed 1234\n";
    file << "emin = 10\n";
    file << "emax = 1000\n";
    file << "esize = 8\n";
    file << "hkpc = 3.5\n";
    file << "discsizepc = 120\n";
    file << "rgkpc = 18\n";
    file << "sunkpc = 8.2\n";
    file << "d0h = 0.5\n";
    file << "e0 = 500\n";
    file << "delta = 0.42\n";
    file << "ddelta = -0.2\n";
    file << "s = 0.15\n";
    file << "eb = 250\n";
    file << "a = 2.1\n";
    file << "b = 4.9\n";
    file << "r1kpc = 0.4\n";
    file << "injslope = 2.2\n";
    file << "injslopesigma = 0.08\n";
    file << "injemax = 50000\n";
    file << "injemaxsigmadex = 0.35\n";
    file << "injemaxmingev = 1000\n";
    file << "injemaxmaxgev = 3000000\n";
    file << "efficiency = 0.18\n";
    file << "bmug = 3.2\n";
    file << "urad = 0.7\n";
    file << "snrateyr = 0.015\n";
    file << "timestep = 0.5\n";
    file << "maxtimemyr = 7\n";
    file << "pid = He\n";
    file << "varyenergy = true\n";
    file << "varyslope = yes\n";
    file << "spiralmodel = Steiman2010\n";
    file << "transportmodel = diffusionlosses\n";
    file << "injectionmodel = galacticrandom\n";
  }

  const core::Input in(path);
  EXPECT_EQ(in.simname(), "parsed_run");
  EXPECT_EQ(in.seed(), 1234u);
  EXPECT_DOUBLE_EQ(in.E_min(), 10.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.E_max(), 1000.0 * cgs::GeV);
  EXPECT_EQ(in.E_size(), 8u);
  EXPECT_DOUBLE_EQ(in.H(), 3.5 * cgs::kpc);
  EXPECT_DOUBLE_EQ(in.h(), 120.0 * cgs::pc);
  EXPECT_DOUBLE_EQ(in.R_g(), 18.0 * cgs::kpc);
  EXPECT_DOUBLE_EQ(in.R_sun(), 8.2 * cgs::kpc);
  EXPECT_DOUBLE_EQ(in.D0_over_H(), 0.5 * cgs::kpc / cgs::Myr);
  EXPECT_DOUBLE_EQ(in.E_0(), 500.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.delta(), 0.42);
  EXPECT_DOUBLE_EQ(in.ddelta(), -0.2);
  EXPECT_DOUBLE_EQ(in.s(), 0.15);
  EXPECT_DOUBLE_EQ(in.E_b(), 250.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.a(), 2.1);
  EXPECT_DOUBLE_EQ(in.b(), 4.9);
  EXPECT_DOUBLE_EQ(in.R_1(), 0.4 * cgs::kpc);
  EXPECT_DOUBLE_EQ(in.injSlope(), 2.2);
  EXPECT_DOUBLE_EQ(in.injSlopeSigma(), 0.08);
  EXPECT_DOUBLE_EQ(in.injEmax(), 50000.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.injEmaxSigmaDex(), 0.35);
  EXPECT_DOUBLE_EQ(in.injEmaxMin(), 1000.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.injEmaxMax(), 3000000.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.injEfficiency(), 0.18);
  EXPECT_DOUBLE_EQ(in.B_field(), 3.2 * cgs::microgauss);
  EXPECT_DOUBLE_EQ(in.U_rad(), 0.7 * cgs::eV / cgs::cm3);
  EXPECT_DOUBLE_EQ(in.sn_rate(), 0.015 / cgs::year);
  EXPECT_DOUBLE_EQ(in.time_step(), 0.5 * cgs::year);
  EXPECT_DOUBLE_EQ(in.max_time(), 7.0 * cgs::Myr);
  EXPECT_EQ(in.pid(), core::He);
  EXPECT_TRUE(in.doVaryEnergy());
  EXPECT_TRUE(in.doVarySlope());
  EXPECT_EQ(in.spiralModel(), SpiralModel::Steiman2010);
  EXPECT_EQ(in.transportModel(), TransportModel::DiffusionLosses);
  EXPECT_EQ(in.injectionModel(), InjectionModel::GalacticRandom);

  EXPECT_EQ(std::remove(path.c_str()), 0);
}

TEST(InputParsing, DerivesSimnameFromFilenameWhenMissing) {
  const std::string path = "/tmp/gryphon_config_derived.txt";
  {
    std::ofstream file(path);
    ASSERT_TRUE(file.is_open());
    file << "emin = 10\n";
    file << "emax = 100\n";
    file << "injemax = 1000\n";
  }

  const core::Input in(path);
  EXPECT_EQ(in.simname(), "gryphon_config_derived");

  EXPECT_EQ(std::remove(path.c_str()), 0);
}

TEST(InputParsing, RejectsUnknownParameter) {
  const auto path = makeTempFilePath();
  {
    std::ofstream file(path);
    ASSERT_TRUE(file.is_open());
    file << "mystery_parameter = 42\n";
  }

  EXPECT_THROW({ core::Input in(path); }, std::invalid_argument);
  EXPECT_EQ(std::remove(path.c_str()), 0);
}

TEST(InputParsing, RepositoryExampleInputFileParses) {
  const core::Input in("data/example_input.txt");

  EXPECT_EQ(in.simname(), "example_run");
  EXPECT_EQ(in.seed(), 12345u);
  EXPECT_EQ(in.E_size(), 64u);
  EXPECT_EQ(in.injectionModel(), InjectionModel::GalacticRandom);
  EXPECT_EQ(in.transportModel(), TransportModel::DiffusionLosses);
}

TEST(InputParsing, ParsesPwnInjectionModel) {
  const auto path = makeTempFilePath();
  {
    std::ofstream file(path);
    ASSERT_TRUE(file.is_open());
    file << "emin = 10\n";
    file << "emax = 100\n";
    file << "injectionmodel = pwn\n";
    file << "p0ms = 120\n";
    file << "sigmap0ms = 35\n";
    file << "pwnrandominitialperiod = false\n";
    file << "pwnalpha1 = 1.4\n";
    file << "pwnalpha2 = 2.8\n";
    file << "pwnebreakgev = 250\n";
    file << "pwnemingev = 2\n";
  }

  const core::Input in(path);
  EXPECT_EQ(in.injectionModel(), InjectionModel::PWN);
  EXPECT_DOUBLE_EQ(in.pwnP0(), 120.0 * cgs::msec);
  EXPECT_DOUBLE_EQ(in.pwnSigmaP0(), 35.0 * cgs::msec);
  EXPECT_FALSE(in.pwnRandomInitialPeriod());
  EXPECT_DOUBLE_EQ(in.pwnAlpha1(), 1.4);
  EXPECT_DOUBLE_EQ(in.pwnAlpha2(), 2.8);
  EXPECT_DOUBLE_EQ(in.pwnEbreak(), 250.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.pwnEmin(), 2.0 * cgs::GeV);

  EXPECT_EQ(std::remove(path.c_str()), 0);
}

TEST(InputParsing, ParsesRandomEmaxInjectionModel) {
  const auto path = makeTempFilePath();
  {
    std::ofstream file(path);
    ASSERT_TRUE(file.is_open());
    file << "emin = 10\n";
    file << "emax = 100\n";
    file << "injectionmodel = randomemax\n";
    file << "injemax = 300000\n";
    file << "injemaxsigmadex = 0.4\n";
    file << "injemaxmingev = 1000\n";
    file << "injemaxmaxgev = 3000000\n";
  }

  const core::Input in(path);
  EXPECT_EQ(in.injectionModel(), InjectionModel::RandomEmax);
  EXPECT_DOUBLE_EQ(in.injEmax(), 300000.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.injEmaxSigmaDex(), 0.4);
  EXPECT_DOUBLE_EQ(in.injEmaxMin(), 1000.0 * cgs::GeV);
  EXPECT_DOUBLE_EQ(in.injEmaxMax(), 3000000.0 * cgs::GeV);

  EXPECT_EQ(std::remove(path.c_str()), 0);
}

TEST(InputParsing, RejectsExplicitPwnCutoffParameter) {
  const auto path = makeTempFilePath();
  {
    std::ofstream file(path);
    ASSERT_TRUE(file.is_open());
    file << "injectionmodel = pwn\n";
    file << "pwnecutgev = 10000\n";
  }

  EXPECT_THROW({ core::Input in(path); }, std::invalid_argument);
  EXPECT_EQ(std::remove(path.c_str()), 0);
}

}  // namespace gryphon
