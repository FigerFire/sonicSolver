#!/usr/bin/env python3
"""Gate B guards for the Program decomposition phase.

Two independent guards, selected with --mode:

  scope       §15/§16/§17/§35 resolvedSystemConsumerScope
              Consumers that only need one compile product must not include the
              ResolvedSimulationSystem aggregate header.  Re-introducing that
              include is how the God object grows back.

  fieldcount  §20/§35 resolvedSystemFieldCount
              The aggregate's declared field count may only decrease.  Raising
              the baseline requires an explicit architectural decision, so this
              guard fails instead of silently accepting new state on the
              aggregate.
"""

import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SYSTEM_DIR = ROOT / "src/solver/system"
AGGREGATE_HEADER = "SF_resolvedSimulationSystem.h"
AGGREGATE = SYSTEM_DIR / AGGREGATE_HEADER

# §20 baseline: number of declared data members of ResolvedSimulationSystem.
# Lower it whenever the aggregate shrinks; never raise it silently.
# §20 baseline: 10 = 原 10 个 compile product。本轮把 solver-family 的
# `formulation` 字段换成两个 typed compile product（CompiledStateRealization
# 与 CouplingReport）：语义 authority 净减少一个 family 开关，但字段数 +1。
FIELD_BASELINE = 11

# Consumers that were narrowed to the specific compile products they need.
# Each entry: file -> documented scope, used in the failure message.
NARROWED_CONSUMERS = {
    "src/solver/run/SF_planExecutor.h": "CompiledSolvePlan",
    "src/solver/run/SF_planExecutor.cpp": "CompiledSolvePlan",
    "src/solver/system/SF_stateRealizer.h": "ExecutableEquationSystem + RuntimeRequirements",
    "src/solver/system/SF_stateRealizer.cpp": "ExecutableEquationSystem + RuntimeRequirements",
    "src/app/application/SF_environment.h": "RuntimeRequirements + mesh/IBM runtime input",
    "src/app/application/SF_environment.cpp": "RuntimeRequirements + mesh/IBM runtime input",
    "src/solver/algorithm/SF_singleFluidStepper.h":
        "ExecutableEquationSystem + CompiledNumericalSystem + CompiledSolvePlan + RuntimeRequirements",
    "src/solver/algorithm/SF_singleFluidStepper.cpp":
        "ExecutableEquationSystem + CompiledNumericalSystem + CompiledSolvePlan + RuntimeRequirements",
    "src/solver/algorithm/eulerian/SF_eulerianStepper.h":
        "ExecutableEquationSystem + CompiledSolvePlan + RuntimeRequirements",
    "src/solver/algorithm/eulerian/SF_eulerianStepper.cpp":
        "ExecutableEquationSystem + CompiledSolvePlan + RuntimeRequirements",
}

# Gate A structure: the native case is decoded straight into typed sections.
# The adapter must own that typed bundle instead of an OpenFOAM-shaped
# pseudo-path document map.
ADAPTER_HEADER = ROOT / "src/app/application/model/compatibility/SF_compatibility.h"
TYPED_SECTION_MEMBERS = [
    "runtime",
    "output",
    "parallel",
    "algorithm",
    "numerics",
    "thermoDynamics",
    "phaseSystem",
    "phaseChange",
    "turbulence",
    "ibm",
    "ilw",
    "gravity",
    "mrf",
    "wallHeat",
]
# Symbols that would mean the OpenFOAM document round-trip came back.
FORBIDDEN_SYMBOLS = [
    "parseOpenFOAM",
    "documents_",
    "isFoamControlDictPath",
    "documentFile(",
    "dictionaryPath(",
]
# Typed section availability helpers the decode stage must be able to query.
REQUIRED_ADAPTER_TOKENS = [
    "struct NativeCaseSections",
    "sections_",
    "const Model::FieldDescriptor* field(",
    "decodeNativeCase",
]

INCLUDE_RE = re.compile(r'^\s*#\s*include\s+"([^"]+)"', re.MULTILINE)


def aggregate_fields():
    """Return declared data member names of ResolvedSimulationSystem."""
    text = AGGREGATE.read_text(errors="replace")
    match = re.search(
        r"struct\s+ResolvedSimulationSystem\s*\{(.*?)\n\};",
        text,
        re.DOTALL,
    )
    if not match:
        sys.exit(f"ERROR scope: cannot locate ResolvedSimulationSystem in {AGGREGATE}")
    body = match.group(1)
    fields = []
    for raw in body.splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        # A declaration line ends with ';' and is not a nested struct/function.
        if not line.endswith(";"):
            continue
        if "(" in line or "struct " in line:
            continue
        fields.append(line.rstrip(";").strip())
    return fields


def check_scope():
    failures = []

    for relative in NARROWED_CONSUMERS:
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"{relative}: narrowed consumer is missing")
            continue
        text = path.read_text(errors="replace")
        if f'"{AGGREGATE_HEADER}"' in text:
            failures.append(
                f"{relative}: narrowed consumer includes the {AGGREGATE_HEADER} "
                f"aggregate again (allowed scope: {NARROWED_CONSUMERS[relative]})"
            )

    # Gate A structure: one typed section per native semantic object, and no
    # OpenFOAM-shaped pseudo-path document map anywhere in the IO path.
    header = ADAPTER_HEADER.read_text(errors="replace")
    for token in REQUIRED_ADAPTER_TOKENS:
        if token not in header:
            failures.append(
                f"{ADAPTER_HEADER.relative_to(ROOT)}: missing typed native IO "
                f"structure '{token}'"
            )
    for section in TYPED_SECTION_MEMBERS:
        if not re.search(rf"\b{section}\b\s*[;{{]", header):
            failures.append(
                f"{ADAPTER_HEADER.relative_to(ROOT)}: typed section slot "
                f"'{section}' is missing"
            )

    for relative, _scope in sorted(NARROWED_CONSUMERS.items()) + sorted(
        (p, "") for p in [
            "src/app/application/model/compatibility/SF_compatibility.h",
            "src/app/application/model/compatibility/SF_compatibility.cpp",
            "src/app/application/model/compatibility/SF_nativeDecode.cpp",
            "src/app/application/model/SF_nativeBinding.cpp",
        ]
    ):
        path = ROOT / relative
        if not path.is_file():
            continue
        text = path.read_text(errors="replace")
        for symbol in FORBIDDEN_SYMBOLS:
            if symbol in text:
                failures.append(
                    f"{relative}: OpenFOAM document round-trip symbol "
                    f"'{symbol}' is back in the native IO path"
                )

    # The execution environment must not infer runtime services from solver
    # identity.
    environment = ROOT / "src/app/application/SF_environment.cpp"
    if environment.is_file() and "solverConfig.numerics.solver" in environment.read_text(
        errors="replace"
    ):
        failures.append(
            "src/app/application/SF_environment.cpp: runtime services inferred "
            "from solver identity"
        )

    return failures


def check_fieldcount():
    fields = aggregate_fields()
    if len(fields) > FIELD_BASELINE:
        return [
            f"{AGGREGATE.relative_to(ROOT)}: field count grew to {len(fields)} "
            f"(baseline {FIELD_BASELINE}). The §20 guard only allows the "
            "aggregate to shrink; new state belongs on a typed compile product.\n"
            "  declared members: " + ", ".join(fields)
        ]
    if len(fields) < FIELD_BASELINE:
        print(
            f"NOTE resolved-system field count decreased to {len(fields)} "
            f"(baseline {FIELD_BASELINE}); lower FIELD_BASELINE."
        )
    return []


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["scope", "fieldcount"], required=True)
    args = parser.parse_args()

    failures = check_scope() if args.mode == "scope" else check_fieldcount()
    if failures:
        for failure in failures:
            print(f"FAIL {args.mode}: {failure}")
        return 1
    if args.mode == "scope":
        print(
            "resolved-system consumers stay inside their compile product scope "
            f"({len(NARROWED_CONSUMERS)} files) and the native IO path has no "
            "OpenFOAM document round-trip"
        )
    else:
        print("resolved-system field count stays within the §20 baseline")
    return 0


if __name__ == "__main__":
    sys.exit(main())
