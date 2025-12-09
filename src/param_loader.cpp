#include "param_loader.hpp"

namespace param_loader {

    template std::string ParamLoader::get_params<std::string>(
    const std::string&, const std::string&) const;
template int ParamLoader::get_params<int>(const std::string&,
                                         const std::string&) const;
template float ParamLoader::get_params<float>(const std::string&,
                                             const std::string&) const;
template double ParamLoader::get_params<double>(const std::string&,
                                             const std::string&) const;
template bool ParamLoader::get_params<bool>(const std::string&,
                                           const std::string&) const;

YAML::Node ParamLoader::navigate_node(const std::string& node_,
                                     const std::string& parameter) const {
  YAML::Node node = config["_CFD_PARAMETERS_"][node_];
  auto param_node = node[parameter];
  if (!param_node) {
    std::cerr << "Parameter path is invalid: " << parameter << std::endl;
    throw std::runtime_error("Invalid parameter path");
  }
  return param_node;
}

template <typename T>
T ParamLoader::get_params(const std::string& node_,
                         const std::string& parameter) const {
  try {
    YAML::Node node = navigate_node(node_, parameter);
    // handle for enum class
    if constexpr (std::is_enum<T>::value) {
      return static_cast<T>(node.as<int>());
    } else {
      return node.as<T>();
    }
  } catch (const YAML::BadConversion& e) {
    std::cerr << "Bad conversion: " << e.what() << std::endl;
    return T();   // Return default-constructed value on failure
  }
}

template <typename T>
std::vector<T> ParamLoader::get_params(const std::string& node_,
                                      const std::string& parameter,
                                      char delimiter) const {
  try {
    YAML::Node node = navigate_node(node_, parameter);
    std::string param_value = node.as<std::string>();
    std::vector<T> values;
    std::istringstream iss(param_value);
    std::string token;
    while (std::getline(iss, token, delimiter)) {
      std::stringstream convert(token);
      T value;
      if (!(convert >> value)) {
        std::cerr << "Conversion failed for token: " << token << std::endl;
        continue;
      }
      values.push_back(value);
    }
    return values;
  } catch (const YAML::BadConversion& e) {
    std::cerr << "Bad conversion: " << e.what() << std::endl;
    return {};   // Return empty vector on failure
  }
}
}   // namespace param_loader
