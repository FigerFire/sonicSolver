#pragma once
/// WHICH occurrence selection is independent of equation math and execution order.
#include <string>
#include <map>
#include <vector>
namespace SF::System {
struct NumericalBinding {
    std::string equation;
    std::string method;
    std::vector<std::string> inputs;
    std::string occurrence;
    /// Selected numerical parameters; never inferred from the equation name.
    std::map<std::string,double> parameters;
};

}
