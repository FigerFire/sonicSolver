set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(QT NAMES Qt5 Qt6 COMPONENTS Widgets QUIET)
if(NOT QT_FOUND)
    message(WARNING "Qt5/Qt6 Widgets was not found; sonicGui will not be built.")
    return()
endif()

find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)

# macOS AMFI 可能拒绝带下载 provenance 的 Homebrew rcc。默认仍使用 Qt 自带
# Qt::rcc；仅在用户显式提供时替换工具路径，便于受限机器完成可复现构建。
set(SONIC_QT_RCC_EXECUTABLE "" CACHE FILEPATH
    "Optional trusted Qt rcc executable used by AUTORCC")
if(SONIC_QT_RCC_EXECUTABLE AND TARGET Qt${QT_VERSION_MAJOR}::rcc)
    if(NOT EXISTS "${SONIC_QT_RCC_EXECUTABLE}")
        message(FATAL_ERROR
            "SONIC_QT_RCC_EXECUTABLE does not exist: ${SONIC_QT_RCC_EXECUTABLE}")
    endif()
    set_target_properties(Qt${QT_VERSION_MAJOR}::rcc PROPERTIES
        IMPORTED_LOCATION "${SONIC_QT_RCC_EXECUTABLE}"
        IMPORTED_LOCATION_RELEASE "${SONIC_QT_RCC_EXECUTABLE}"
        IMPORTED_LOCATION_RELWITHDEBINFO "${SONIC_QT_RCC_EXECUTABLE}"
        IMPORTED_LOCATION_MINSIZEREL "${SONIC_QT_RCC_EXECUTABLE}"
        IMPORTED_LOCATION_DEBUG "${SONIC_QT_RCC_EXECUTABLE}"
    )
endif()

find_package(VTK QUIET COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersGeometry
    GUISupportQt
    InteractionStyle
    IOGeometry
    IOXML
    RenderingAnnotation
    RenderingCore
    RenderingOpenGL2
)

if(VTK_FOUND)
    set(SONIC_GUI_VTK_FOUND ON)
    message(STATUS "sonicGui: VTK preview enabled")
else()
    set(SONIC_GUI_VTK_FOUND OFF)
    message(STATUS "sonicGui: VTK not found; building without 3D preview")
endif()

add_subdirectory(IO)
add_subdirectory(Geometry)
add_subdirectory(ViewerControls)
add_subdirectory(VTKView)
add_subdirectory(Resources)
add_subdirectory(hotKey)
add_subdirectory(GUI)
add_subdirectory(Controller)
add_subdirectory(Splash)

add_library(SF_gui INTERFACE)
target_sources(SF_gui INTERFACE SF_gui.h)
target_include_directories(SF_gui INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(SF_gui INTERFACE
    SF_guiController
)

add_subdirectory(Application)
