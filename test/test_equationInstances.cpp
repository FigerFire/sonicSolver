#include "app/application/model/SF_model.h"

#include <iostream>

int main() {
    SF::Model::Description model;
    model.objects.push_back({
        "liquid","solver","equationSystem",{},
        {{"preset","eulerianEulerian"}}, {}});
    model.objects.push_back({
        "gas","solver","equationSystem",{},
        {{"preset","eulerianEulerian"}}, {}});
    const auto instances =
        SF::Application::ModelLoader::equationInstances(model);
    if (instances.size() != 2
        || instances[0].name != "liquid" || instances[1].name != "gas"
        || instances[0].type != "eulerian"
        || instances[1].type != "eulerian") {
        std::cerr << "Independent same-type equation instances were not preserved.\n";
        return 1;
    }
    model.objects[1].name = "liquid";
    try {
        (void)SF::Application::ModelLoader::equationInstances(model);
    } catch (const std::runtime_error&) {
        return 0;
    }
    std::cerr << "Duplicate equation-system instance name was accepted.\n";
    return 1;
}
