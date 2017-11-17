atlas_subdir( TopoStudies )

set( extra_deps )
set( extra_libs )

atlas_depends_on_subdirs(
    PUBLIC
    Event/xAOD/xAODBase
    Control/xAODRootAccess
    Control/AthToolSupport/AsgTools
    Event/xAOD/xAODEventInfo
    Event/xAOD/xAODTruth
    Trigger/TrigConfiguration/TrigConfxAOD
    Trigger/TrigAnalysis/TrigDecisionTool
    Event/xAOD/xAODTrigMissingET
    Event/xAOD/xAODMissingET
    Event/xAOD/xAODJet
    Event/xAOD/xAODTrigger
    ${extra_deps}
)

add_definitions(-DR21ANA)

find_package( ROOT COMPONENTS Gpad Graf Graf3d Core Tree MathCore Hist Tree RIO )

atlas_add_library( TopoStudiesLib
    TopoStudies/*.h Root/*.cxx
    PUBLIC_HEADERS TopoStudies
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    LINK_LIBRARIES xAODBase xAODRootAccess xAODEventInfo TrigConfxAODLib
    xAODTrigger AsgTools xAODJet TrigDecisionToolLib ${extra_libs}
)

atlas_add_executable( runTopoStudies "util/runTopoStudies.cxx"
    INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    LINK_LIBRARIES ${ROOT_LIBRARIES} xAODRootAccess xAODEventInfo
    xAODCore xAODBase TopoStudiesLib ${extra_libs}
)
