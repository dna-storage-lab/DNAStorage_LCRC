cmake_minimum_required(VERSION 3.10)

project(CompositeCode)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${PROJECT_SOURCE_DIR}/src)

link_directories(lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)


add_executable(IdentificationCons_OligoPool
    src/identification_consensus_oligoPool.cpp
    src/edlib.cpp
)
target_link_libraries(IdentificationCons_OligoPool lcrc pthread)


add_executable(ErrorCorrection_OligoPool
    src/error_correction_oligoPool.cpp
)
target_link_libraries(ErrorCorrection_OligoPool dec pthread)

add_executable(Decoder_RealTime
    src/Decoder_realTime.cpp
    src/edlib.cpp
    src/thpool.cpp
)
target_link_libraries(Decoder_RealTime lcrc dec pthread)

