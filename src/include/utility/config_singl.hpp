#ifndef CONFIG_SINGLTONE_HPP
#define CONFIG_SINGLTONE_HPP

//#include "verify.hpp"

#include <string>

namespace config_common {

// config singleton-wrap
template<class Config>
struct config {
  static std::string dirnameOf(std::string const & fname) {
    size_t pos = fname.find_last_of("\\/");
    return (std::string::npos == pos) ? "" : fname.substr(0, pos);
  }

  static void create_instance(std::string const & filename) {
    load(inner_cfg(), filename);
    is_initialized() = true;
  }

  static Config const& get() {
    //VERIFY_MSG(is_initialized(), "Config not initialized");
    return inner_cfg();
  }

  static Config& get_writable() {
    //VERIFY_MSG(is_initialized(), "Config not initialized");
    return inner_cfg();
  }

private:
  static Config& inner_cfg() {
    static Config config;
    return config;
  }

  static bool& is_initialized() {
    static bool is_initialized = false;
    return is_initialized;
  }
};

}

#endif