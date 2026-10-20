#pragma once

/// @file SF_interfaces.h
/// @brief External compatibility umbrella; production code includes exact ports.
///
/// Internal production includes are rejected by the architecture checker. This
/// header remains temporarily for downstream source compatibility only.

#include "SF_boundaryPipeline.h"
#include "SF_parallelCoordinatorContract.h"
#include "SF_solverStepper.h"
