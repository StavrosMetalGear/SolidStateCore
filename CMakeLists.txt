cmake_minimum_required(VERSION 3.15)
project(SolidStateCore LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Dependencies via FetchContent ──────────────────────────────────────────────
include(FetchContent)

FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
)
set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.91.8
)
FetchContent_MakeAvailable(imgui)

FetchContent_Declare(
    implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG        v0.16
)
FetchContent_MakeAvailable(implot)

# ── ImGui library (OpenGL3 + GLFW backend) ────────────────────────────────────
add_library(imgui_impl STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
)
target_include_directories(imgui_impl PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
    ${implot_SOURCE_DIR}
)
target_link_libraries(imgui_impl PUBLIC glfw)

# ── Find OpenGL ────────────────────────────────────────────────────────────────
find_package(OpenGL REQUIRED)

# ── Static library: SolidStatePhysics (no GUI) ────────────────────────────────
add_library(SolidStatePhysics STATIC
    SolidStateSystem.cpp
    NumericalSolverSS.cpp
    solidstate_c_api.cpp
)
target_include_directories(SolidStatePhysics PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Eigen
)
if(MSVC)
    target_compile_definitions(SolidStatePhysics PUBLIC _USE_MATH_DEFINES NOMINMAX)
endif()

# ── GUI executable: SolidStateCore ─────────────────────────────────────────────
add_executable(SolidStateCore
    main.cpp
    GuiApp.cpp
)
target_link_libraries(SolidStateCore PRIVATE
    SolidStatePhysics
    imgui_impl
    OpenGL::GL
)
target_include_directories(SolidStateCore PRIVATE
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
    ${implot_SOURCE_DIR}
)
if(MSVC)
    target_compile_definitions(SolidStateCore PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_precompile_headers(SolidStateCore PRIVATE pch.h)
endif()
