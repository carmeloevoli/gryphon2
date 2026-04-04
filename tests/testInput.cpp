#include <string>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {

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

TEST(InputValidation, RejectsInvalidInjectionCutoff) {
  core::Input in;
  in.set_injectionModel(InjectionModel::SinglePowerLaw);
  in.set_refEnergy(1.0 * cgs::TeV);
  in.set_injEmax(0.5 * cgs::TeV);
  EXPECT_THROW(in.validate(), std::invalid_argument);
}

}  // namespace gryphon
