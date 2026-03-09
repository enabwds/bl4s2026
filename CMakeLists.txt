cmake_minimum_required(VERSION 3.16)
project(BL4S_Sim)

# ── Find GEANT4 ───────────────────────────────────────────────────────────────
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# ── Source files ──────────────────────────────────────────────────────────────
set(SOURCES
    src/DetectorConstruction.cc
    src/DetectorMessenger.cc
    src/CalorimeterSD.cc
    src/ActionInitialisation.cc
)

# ── Build executable ──────────────────────────────────────────────────────────
add_executable(BL4S_sim BL4S_sim.cc ${SOURCES})

target_include_directories(BL4S_sim PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIRS}          # explicit – ensures all .cc files see G4 headers
)

target_link_libraries(BL4S_sim ${Geant4_LIBRARIES})

# ── Copy macros to build directory ───────────────────────────────────────────
file(GLOB MACROS ${PROJECT_SOURCE_DIR}/macros/*.mac)
foreach(_mac ${MACROS})
    get_filename_component(_mac_name ${_mac} NAME)
    configure_file(${_mac} ${PROJECT_BINARY_DIR}/macros/${_mac_name} COPYONLY)
endforeach()
