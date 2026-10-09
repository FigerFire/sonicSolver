cmake_minimum_required(VERSION 3.15)
project(sonicFluid LANGUAGES C CXX)

# ============================================================
#  1. 编译器设置
# ============================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    execute_process(COMMAND xcrun --show-sdk-path
                    OUTPUT_VARIABLE APPLE_SDK_PATH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(CMAKE_OSX_SYSROOT ${APPLE_SDK_PATH})
    # 物理状态诊断依赖 IEEE-754 NaN/Inf 语义；-ffast-math 会删除
    # std::isfinite/fail-fast 分支，因此求解器全局禁止该选项。
    add_compile_options(-O3)
    message(STATUS "Detected macOS: SDK = ${APPLE_SDK_PATH}")
endif()

if(UNIX AND NOT APPLE)
    add_compile_options(-O3 -march=native -fopenmp)
    message(STATUS "Detected Linux")
endif()

# ============================================================
#  2. 公共头文件 + 辅助宏
# ============================================================
add_library(SF_headers INTERFACE)
target_include_directories(SF_headers INTERFACE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/
)

macro(set_lib_output _target _subdir)
    set_target_properties(${_target} PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/${_subdir}"
    )
endmacro()

# ============================================================
#  3. 源码域 — 模块 target 仍由各自 CMakeLists.txt 独立定义
# ============================================================
add_subdirectory(src/core)
add_subdirectory(src/models)
add_subdirectory(src/infrastructure)
add_subdirectory(src/methods)
add_subdirectory(src/solver)

# ============================================================
#  4. OpenMP (可选)
# ============================================================
find_package(OpenMP)

# ============================================================
#  5. 前端入口 (可选)
# ============================================================
option(BUILD_CLI "Build the sonicSolver command-line application" ON)
option(BUILD_GUI "Build the sonicGui desktop application" ON)

add_subdirectory(src/app)

if(NOT OpenMP_CXX_FOUND)
    message(WARNING "OpenMP not found — single core only.")
endif()

# ============================================================
#  6. 测试
# ============================================================
option(BUILD_TESTS "Build unit test targets" OFF)

if(BUILD_TESTS)
    enable_testing()
    find_package(MPI QUIET COMPONENTS CXX)
    if(EXISTS "${CMAKE_SOURCE_DIR}/test/test_registryIO.cpp")
        add_executable(test_registryIO test/test_registryIO.cpp)
        target_link_libraries(test_registryIO PRIVATE SF_io SF_resultWriter)
        target_compile_definitions(test_registryIO PRIVATE SF_TEST_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
        add_test(NAME registryIO COMMAND test_registryIO)
    endif()

    add_executable(test_numericalFlux test/test_numericalFlux.cpp)
    target_link_libraries(test_numericalFlux
        PRIVATE SF_flux SF_weno3 SF_weno5 SF_teno5 SF_weno7 SF_physics)
    add_test(NAME numericalFluxContract COMMAND test_numericalFlux)

    add_executable(test_pisoArchitecture test/test_pisoArchitecture.cpp)
    target_link_libraries(test_pisoArchitecture PRIVATE SF_solverSystem SF_physics)
    add_test(NAME pisoArchitecture COMMAND test_pisoArchitecture)
    if(MPI_CXX_FOUND AND TARGET SF_mpi)
        add_executable(test_formulaExecution test/test_formulaExecution.cpp)
        target_link_libraries(test_formulaExecution PRIVATE
            SF_solverSystem SF_run SF_linearAlgebra SF_mpi MPI::MPI_CXX)
        add_test(NAME formulaExecutionContract COMMAND test_formulaExecution)
    endif()
    add_executable(test_rhieChow test/test_rhieChow.cpp)
    target_link_libraries(test_rhieChow PRIVATE SF_headers)
    add_test(NAME rhieChowFaceCoupling COMMAND test_rhieChow)

    add_executable(test_fixedTimePressure test/test_fixedTimePressure.cpp)
    target_link_libraries(test_fixedTimePressure PRIVATE SF_headers)
    add_test(NAME fixedTimePressureContract COMMAND test_fixedTimePressure)

    add_executable(test_runtimeBinding test/test_runtimeBinding.cpp)
    add_executable(test_standaloneScalar test/test_standaloneScalar.cpp)
    target_link_libraries(test_standaloneScalar PRIVATE SF_application)
    add_test(NAME standaloneScalarExecution COMMAND test_standaloneScalar)
    add_executable(test_multiScalar test/test_multiScalar.cpp)
    target_link_libraries(test_multiScalar PRIVATE SF_application)
    add_test(NAME multiScalarExecution COMMAND test_multiScalar)
    add_executable(test_coupledScalar test/test_coupledScalar.cpp)
    target_link_libraries(test_coupledScalar PRIVATE SF_application)
    add_test(NAME coupledScalarExecution COMMAND test_coupledScalar)
    add_executable(test_fluidScalar test/test_fluidScalar.cpp)
    target_link_libraries(test_fluidScalar PRIVATE SF_application)
    add_test(NAME fluidScalarExecution COMMAND test_fluidScalar)
    target_link_libraries(test_runtimeBinding PRIVATE SF_solverAlgorithm SF_solverSystem)
    add_test(NAME runtimeBindingLifetime COMMAND test_runtimeBinding)

    add_executable(test_explicitStages test/test_explicitStages.cpp)
    target_link_libraries(test_explicitStages PRIVATE SF_solverAlgorithm SF_solverSystem)
    add_test(NAME explicitStageMathematics COMMAND test_explicitStages)

    add_executable(test_mixedTemporalComposition test/test_mixedTemporalComposition.cpp)
    target_link_libraries(test_mixedTemporalComposition PRIVATE SF_solverSystem SF_run)
    add_test(NAME mixedTemporalComposition COMMAND test_mixedTemporalComposition)
    add_executable(test_nativeEulerian test/test_nativeEulerian.cpp)
    target_link_libraries(test_nativeEulerian PRIVATE SF_solverSystem SF_physics SF_turbulence SF_run)
    add_test(NAME nativeEulerianExecution COMMAND test_nativeEulerian)
    add_executable(test_eulerianAssembly test/test_eulerianAssembly.cpp)
    target_link_libraries(test_eulerianAssembly PRIVATE SF_solverSystem)
    add_test(NAME eulerianAssemblyAuthority COMMAND test_eulerianAssembly)
    add_executable(test_eulerianTurbulenceKernel test/test_eulerianTurbulenceKernel.cpp)
    target_link_libraries(test_eulerianTurbulenceKernel PRIVATE SF_solverAlgorithm)
    add_test(NAME eulerianTurbulenceMatrixKernel COMMAND test_eulerianTurbulenceKernel)
    add_executable(test_eulerianTurbulence test/test_eulerianTurbulence.cpp)
    target_link_libraries(test_eulerianTurbulence PRIVATE SF_solverSystem SF_turbulence SF_run)
    add_test(NAME nativeEulerianTurbulenceContract COMMAND test_eulerianTurbulence)
    add_executable(test_nativeLevelSet test/test_nativeLevelSet.cpp)
    target_link_libraries(test_nativeLevelSet PRIVATE SF_solverSystem SF_physics SF_run)
    add_test(NAME nativeLevelSetExecution COMMAND test_nativeLevelSet)
    add_executable(test_nativeMultiphase test/test_nativeMultiphase.cpp)
    target_link_libraries(test_nativeMultiphase PRIVATE SF_solverSystem SF_physics)
    add_test(NAME nativeMultiphaseCapability COMMAND test_nativeMultiphase)
    add_executable(test_nativeTurbulence test/test_nativeTurbulence.cpp)
    target_link_libraries(test_nativeTurbulence PRIVATE SF_solverSystem SF_turbulence SF_run)
    add_test(NAME nativeTurbulenceExecution COMMAND test_nativeTurbulence)
    add_executable(test_methodObjects test/test_methodObjects.cpp)
    target_link_libraries(test_methodObjects PRIVATE SF_solverSystem SF_run)
    add_test(NAME methodObjectCompilation COMMAND test_methodObjects)

    add_executable(test_providerCatalog test/test_providerCatalog.cpp)
    target_link_libraries(test_providerCatalog PRIVATE SF_solverSystem SF_physics)
    add_test(NAME providerCatalogResolution COMMAND test_providerCatalog)

    add_executable(test_termRecipes test/test_termRecipes.cpp)
    target_link_libraries(test_termRecipes PRIVATE SF_solverSystem SF_physics)
    add_test(NAME termRecipeAuthority COMMAND test_termRecipes)

    add_executable(test_equationInstances test/test_equationInstances.cpp)
    target_link_libraries(test_equationInstances PRIVATE SF_caseAdapter)
    add_test(NAME equationInstanceComposition COMMAND test_equationInstances)

    add_executable(test_nativeIoNoFoamRoundTrip
        test/test_nativeIoNoFoamRoundTrip.cpp)
    target_link_libraries(test_nativeIoNoFoamRoundTrip PRIVATE SF_caseAdapter)
    target_compile_definitions(test_nativeIoNoFoamRoundTrip
        PRIVATE SF_TEST_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    add_test(NAME nativeIoNoFoamRoundTrip COMMAND test_nativeIoNoFoamRoundTrip)

    add_executable(test_programDecomposition test/test_programDecomposition.cpp)
    target_link_libraries(test_programDecomposition PRIVATE SF_solverSystem)
    add_test(NAME programDoesNotOwnRuntimeState COMMAND test_programDecomposition)

    if(EXISTS "${CMAKE_SOURCE_DIR}/test/test_formulationArchitecture.cpp")
        # coupling preset / state realization / compiled HOW 的架构契约。
        add_executable(test_formulationArchitecture
            test/test_formulationArchitecture.cpp)
        target_link_libraries(test_formulationArchitecture
            PRIVATE SF_solverSystem)
        add_test(NAME formulationArchitecture
            COMMAND test_formulationArchitecture)
    endif()

    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(TARGET sonicSolver AND Python3_Interpreter_FOUND)
        add_test(NAME nativeInterfaceNumericalIO
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/check_native_interface.py
                --solver $<TARGET_FILE:sonicSolver> --root ${CMAKE_SOURCE_DIR})
        foreach(model kOmegaSST kEpsilon)
            add_test(NAME eulerian${model}NumericalRegression
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/check_eulerian_turbulence.py
                    --solver $<TARGET_FILE:sonicSolver>
                    --case ${CMAKE_SOURCE_DIR}/test/turbulence/eulerian${model})
        endforeach()
        add_test(NAME eulerianNativeNumericalRegression
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/check_eulerian_native.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/eulerianEulerianCase)
        add_test(NAME singleFluidTurbulenceNumericalRegression
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_single_fluid_turbulence.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/turbulence/singleFluidSST)
        add_test(NAME explicitSolutionStateSelection
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_solution_state_selection.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/pressure/constantDensityPiso)
        add_test(NAME cliGeneratedRecipe
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_cli_init_recipe.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/Sod/sodCase_weno7)
        add_test(NAME sodRusanovNumericalRegression
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_sod_rusanov_regression.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/Sod/sodCase_weno7_t0p2)
        add_test(NAME pisoNumericalRegression
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_piso_regression.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/pressureConstraintPiso)
        add_test(NAME constantDensityPressureRegression
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_constant_density_pressure.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/pressure/constantDensityPiso)
        add_test(NAME constantDensityPoiseuilleRegression
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_poiseuille_pressure.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/pressure/constantDensityPoiseuille)
        add_test(NAME fixedTimePressureRegression
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_fixed_time_pressure.py
                --solver $<TARGET_FILE:sonicSolver>
                --cavity ${CMAKE_SOURCE_DIR}/test/pressure/constantDensityPiso
                --poiseuille ${CMAKE_SOURCE_DIR}/test/pressure/constantDensityPoiseuille)
        add_test(NAME pressureTermProviders
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_pressure_term_providers.py
                --solver $<TARGET_FILE:sonicSolver>
                --case ${CMAKE_SOURCE_DIR}/test/pressure/constantDensityPiso)
    endif()

    if(Python3_Interpreter_FOUND)
        add_test(NAME resolvedSystemConsumerScope
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_program_decomposition.py
                --mode scope)
        add_test(NAME resolvedSystemFieldCount
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tools/check_program_decomposition.py
                --mode fieldcount)
    endif()


    add_executable(test_executionContract test/test_executionContract.cpp)
    target_link_libraries(test_executionContract PRIVATE SF_solverSystem SF_physics)
    add_test(NAME scopeAndStateContracts COMMAND test_executionContract)

    add_executable(test_viscousImmersedWall test/test_viscousImmersedWall.cpp)
    target_link_libraries(test_viscousImmersedWall PRIVATE SF_solverSystem SF_ibm SF_turbulence SF_run)
    add_test(NAME viscousImmersedWall COMMAND test_viscousImmersedWall)
    if(Python3_Interpreter_FOUND AND BUILD_CLI)
        add_test(NAME viscousSSTNativeSchedule
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/check_viscous_immersed_wall.py
                --solver $<TARGET_FILE:sonicSolver> --root ${CMAKE_SOURCE_DIR})
    endif()

    add_executable(test_nativeImmersed test/test_nativeImmersed.cpp)
    target_link_libraries(test_nativeImmersed PRIVATE SF_solverSystem SF_ibm SF_turbulence SF_run)
    add_test(NAME nativeImmersedContributions COMMAND test_nativeImmersed)

    if(EXISTS "${CMAKE_SOURCE_DIR}/test/test_field.cpp")
        add_executable(test_field test/test_field.cpp)
        target_link_libraries(test_field PRIVATE SF_field)
    endif()

    if(EXISTS "${CMAKE_SOURCE_DIR}/test/test_fieldVectorLayout.cpp")
        # 保护 Field 变量组布局 contract：Vector3 的 offset 必须来自 v，
        # 不能假定 RU/RV/RW == 1/2/3。
        add_executable(test_fieldVectorLayout test/test_fieldVectorLayout.cpp)
        target_link_libraries(test_fieldVectorLayout PRIVATE SF_field)
        add_test(NAME fieldVectorLayout
            COMMAND test_fieldVectorLayout)
    endif()

    if(EXISTS "${CMAKE_SOURCE_DIR}/test/test_meshGen.cpp")
        add_executable(test_meshGen test/test_meshGen.cpp)
        target_link_libraries(test_meshGen PRIVATE SF_io SF_mesh)
    endif()

    if(EXISTS "${CMAKE_SOURCE_DIR}/test/test_canonicalFace.cpp")
        add_executable(test_canonicalFace test/test_canonicalFace.cpp)
        target_link_libraries(test_canonicalFace PRIVATE SF_field)
        add_test(NAME canonicalFaceContract COMMAND test_canonicalFace)
    endif()

    add_executable(test_nodalQuadrature test/test_nodalQuadrature.cpp)
    target_link_libraries(test_nodalQuadrature PRIVATE SF_field)
    add_test(NAME nodalQuadrature COMMAND test_nodalQuadrature)
    add_executable(test_surfaceMoments test/test_surfaceMoments.cpp)
    target_link_libraries(test_surfaceMoments PRIVATE SF_ibm)
    add_test(NAME surfaceFirstMoments COMMAND test_surfaceMoments)
    if(TARGET SF_mpi AND MPI_CXX_FOUND)
        add_executable(test_distributedSurfaceMoments test/test_distributedSurfaceMoments.cpp)
        target_link_libraries(test_distributedSurfaceMoments PRIVATE SF_ibm SF_mpi SF_executionRuntime MPI::MPI_CXX)
        add_test(NAME distributedSurfaceFirstMoments2 COMMAND mpirun -np 2 $<TARGET_FILE:test_distributedSurfaceMoments>)
        add_test(NAME distributedSurfaceFirstMoments4 COMMAND mpirun -np 4 $<TARGET_FILE:test_distributedSurfaceMoments>)
    endif()

    add_executable(test_physicalILWPartition test/test_physicalILWPartition.cpp)
    target_link_libraries(test_physicalILWPartition PRIVATE SF_boundary)
    add_test(NAME physicalILWPartition COMMAND test_physicalILWPartition)

    if(EXISTS "${CMAKE_SOURCE_DIR}/test/test_surfaceSchurProjection.cpp")
        add_executable(test_surfaceSchurProjection
            test/test_surfaceSchurProjection.cpp)
        target_link_libraries(test_surfaceSchurProjection PRIVATE SF_ibm)
        add_test(NAME surfaceSchurProjection
            COMMAND test_surfaceSchurProjection)
    endif()

    if(EXISTS "${CMAKE_SOURCE_DIR}/test/test_distributedConstraintLayout.cpp")
        add_executable(test_distributedConstraintLayout
            test/test_distributedConstraintLayout.cpp)
        target_link_libraries(test_distributedConstraintLayout
            PRIVATE SF_linearAlgebra SF_executionRuntime SF_mpi MPI::MPI_CXX)
        add_test(NAME distributedConstraintLayout
            COMMAND mpirun -np 4 $<TARGET_FILE:test_distributedConstraintLayout>)
    endif()

    if(TARGET SF_mpi AND MPI_CXX_FOUND)
        add_executable(test_haloFailureContract test/test_haloFailureContract.cpp)
        target_link_libraries(test_haloFailureContract PRIVATE SF_mpi SF_mesh MPI::MPI_CXX)
        add_test(NAME haloFailureContract
            COMMAND mpirun -np 2 $<TARGET_FILE:test_haloFailureContract>)
        set_tests_properties(haloFailureContract PROPERTIES TIMEOUT 30)
        add_executable(test_distributedPressure test/test_distributedPressure.cpp)
        target_link_libraries(test_distributedPressure PRIVATE SF_linearAlgebra SF_mpi SF_executionRuntime MPI::MPI_CXX)
        add_test(NAME distributedPressurePrimitives
            COMMAND mpirun -np 2 $<TARGET_FILE:test_distributedPressure>)
        set_tests_properties(distributedPressurePrimitives PROPERTIES TIMEOUT 30)
        add_test(NAME distributedPressureFailStop
            COMMAND mpirun -np 2 $<TARGET_FILE:test_distributedPressure> --fail-on-one-rank)
        set_tests_properties(distributedPressureFailStop PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
        if(TARGET sonicSolver AND Python3_Interpreter_FOUND)

            add_test(NAME distributedPressureRegression
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/check_distributed_pressure.py
                    --solver $<TARGET_FILE:sonicSolver>
                    --cavity ${CMAKE_SOURCE_DIR}/test/pressure/constantDensityPiso
                    --poiseuille ${CMAKE_SOURCE_DIR}/test/pressure/constantDensityPoiseuille)
            set_tests_properties(distributedPressureRegression PROPERTIES TIMEOUT 300)
        endif()
    endif()

    if(TARGET SF_mpi AND MPI_CXX_FOUND)
        add_executable(test_peskinPartialRows test/test_peskinPartialRows.cpp)
        target_link_libraries(test_peskinPartialRows PRIVATE SF_ibm SF_executionRuntime SF_mpi MPI::MPI_CXX)
        add_test(NAME peskinPartialRows COMMAND mpirun -np 4 $<TARGET_FILE:test_peskinPartialRows>
            ${CMAKE_SOURCE_DIR}/test/IBM/cylinderFlowPeskin/mesh/geometry/cylinderSurface.stl)
    endif()

    if(TARGET SF_mpi AND MPI_CXX_FOUND
       AND EXISTS "${CMAKE_SOURCE_DIR}/test/test_sparseCanonicalCopy.cpp"
       AND EXISTS "${CMAKE_SOURCE_DIR}/test/test_distributedKKT.cpp"
       AND EXISTS "${CMAKE_SOURCE_DIR}/test/test_distributedSurfaceSchur.cpp")
        add_executable(test_sparseCanonicalCopy
            test/test_sparseCanonicalCopy.cpp)
        target_link_libraries(test_sparseCanonicalCopy
            PRIVATE SF_mpi SF_executionRuntime MPI::MPI_CXX)
        add_test(NAME sparseCanonicalCopy
            COMMAND mpirun -np 4 $<TARGET_FILE:test_sparseCanonicalCopy>)

        add_executable(test_distributedKKT
            test/test_distributedKKT.cpp)
        target_link_libraries(test_distributedKKT
            PRIVATE SF_linearAlgebra SF_mpi SF_executionRuntime MPI::MPI_CXX)
        add_test(NAME distributedMinimalKKT
            COMMAND mpirun -np 2 $<TARGET_FILE:test_distributedKKT>)

        add_executable(test_distributedSurfaceSchur
            test/test_distributedSurfaceSchur.cpp)
        target_link_libraries(test_distributedSurfaceSchur
            PRIVATE SF_ibm SF_mpi SF_executionRuntime MPI::MPI_CXX)
        add_test(NAME distributedSurfaceSchur
            COMMAND mpirun -np 4 $<TARGET_FILE:test_distributedSurfaceSchur>)
    endif()

endif()
