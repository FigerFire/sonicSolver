#pragma once
/// @brief 持有命名标量的唯一物理数组及边界端口；不持有 clock、stage workspace 或执行循环。
#include "core/field/SF_scalarField.h"
#include "core/interfaces/SF_scalarEquationData.h"
#include "core/system/SF_stateViews.h"
#include <memory>
#include <vector>
namespace SF { struct CaseConfig; class Field; namespace State { struct StateBundle; }
namespace System { struct ResolvedSimulationSystem; }
namespace Application::Execution::Detail {
struct ScalarResources {
    std::vector<std::unique_ptr<ScalarField>> arrays;
    std::vector<std::shared_ptr<const FDM::ScalarEquationData>> ports;
    std::vector<FDM::ScalarInstanceData> instances;
    std::vector<System::ResourceView> views;
};
/// @brief 消费 compiled physicalResources；同一 binder 服务 flow+scalar 与 geometry-only。
ScalarResources bindScalarResources(Field&,State::StateBundle&,
    const System::ResolvedSimulationSystem&,const CaseConfig&,int block=0,
    const std::vector<System::ResourceView>& alreadyBound={});
}}
