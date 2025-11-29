cmake_minimum_required(VERSION 3.10)

project(CompositeCode)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${PROJECT_SOURCE_DIR}/src)

link_directories(lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_executable(Gen_LCRC
    src/Generate_LCRC_refer.cpp
)
target_link_libraries(Gen_LCRC lcrc)

add_executable(Decoder_OligoPool
    src/Decoder_oligoPool_300K.cpp
    src/edlib.cpp
    src/thpool.cpp    
)
target_link_libraries(Decoder_OligoPool lcrc dec pthread)

add_executable(Decoder_RealTime
    src/Decoder_oligoPool_realTime_300K.cpp
    src/edlib.cpp
    src/thpool.cpp
)
target_link_libraries(Decoder_RealTime lcrc dec pthread)

