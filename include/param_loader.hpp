#pragma once 

#include <yaml-cpp/yaml.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace param_loader {
class ParamLoader {
 private:
  static bool isInitialized_;
  YAML::Node config;

  /**
   * @brief Constructor to load the YAML configuration file
   *
   */
  explicit ParamLoader(const std::string& filename) {
    try {
      config = YAML::LoadFile(filename);
    } catch (const YAML::BadFile& e) {
      std::cerr << "Failed to open file: " << e.what() << std::endl;
      throw std::runtime_error("Failed to load YAML configuration.");
    }
  }

  /**
   * @brief Navigate to the node in the YAML configuration file
   *
   * @param node
   * @param parameter
   * @return YAML::Node
   */
  YAML::Node navigate_node(const std::string& node,
                          const std::string& parameter) const;

 public:
  /**
   * @brief Deleted copy constructor
   *
   */
  ParamLoader(const ParamLoader&) = delete;

  /**
   * @brief Deleted assignment operator
   *
   * @return ParamLoader&
   */
  ParamLoader& operator=(const ParamLoader&) = delete;

  /**
   * @brief Get the Instance object
   *
   * @param filename
   * @return ParamLoader&
   */
  static ParamLoader& getInstance(const std::string& filename) {
    static ParamLoader instance(filename);
    return instance;
  }

  // Static method to get the singleton instance as a shared pointer
  static std::shared_ptr<ParamLoader>
  getSharedInstance(const std::string& filename) {
    static std::shared_ptr<ParamLoader> instance(new ParamLoader(filename));
    return instance;
  }

  /**
   * @brief Getter method to get parameters from the YAML configuration
   *
   * @tparam T
   * @param node
   * @param parameter
   * @return T
   */
  template <typename T>
  T get_params(const std::string& node, const std::string& parameter) const;

  /**
   * @brief Getter method to get parameters as a vector from the YAML
   * configuration
   *
   * @param node
   * @param parameter
   * @param delimiter
   * @return * template <typename T>
   */
  template <typename T>
  std::vector<T> get_params(const std::string& node,
                           const std::string& parameter, char delimiter) const;

  // operator for std::cout
  friend std::ostream& operator<<(std::ostream& os, const ParamLoader& param) {
    os << param.config;
    return os;
  }
};
}