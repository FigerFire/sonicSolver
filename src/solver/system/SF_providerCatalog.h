#pragma once

#include "SF_runtimeRequirements.h"
#include <string_view>

namespace SF::System {

struct ProviderDescriptor {
    std::string id;
    std::vector<OperationCapability> capabilities;
};

/// Validate the implementation selected by WHICH; never nominate alternatives.
class ProviderCatalog {
public:
    void add(ProviderDescriptor descriptor);
    ResolvedOperationBinding resolve(const ExecutableOperation& operation,
                                     std::string_view selectedProvider) const;
    static ProviderCatalog builtIn();
private:
    std::vector<ProviderDescriptor> descriptors_;
};

} // namespace SF::System
