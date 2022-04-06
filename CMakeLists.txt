cmake_minimum_required(VERSION 3.15)

project(hyperborean-sun)

include(thirdparty/cmake-find-or-download/find-or-download.cmake)

# Download glslang-validator and shaderboi for testing.
find_or_download_if_not_present(SHADERBOI 
    sboi.exe
    "https://github.com/Boscop/ShaderBoi/releases/download/v1.0.0/ShaderBoi_v1.0.0-windows.zip" 
    ShaderBoi/
)
find_or_download_if_not_present(GLSLANG_VALIDATOR 
    glslangValidator.exe
    "https://github.com/KhronosGroup/glslang/releases/download/master-tot/glslang-master-windows-x64-Release.zip" 
    bin/
)

# Obtain glext.h and khrplatform.h.
find_or_download_if_not_present(GLEXT
    glext.h "https://www.khronos.org/registry/OpenGL/api/GL/glext.h"
    ./
)
find_or_download_if_not_present(KHRPLATFORM
    khrplatform.h
    "https://www.khronos.org/registry/EGL/api/KHR/khrplatform.h"
    ./
)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/downloads/KHR")
file(COPY "${CMAKE_CURRENT_BINARY_DIR}/downloads/khrplatform.h" 
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/downloads/KHR")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/downloads")

# Find the shadertoy-to-c-exporter binary.
find_program(SHADERTOY_TO_C NAMES shadertoy-to-c.exe PATHS "${PROJECT_SOURCE_DIR}/bin")

# And we need OpenGL.
find_package(OpenGL REQUIRED)

add_subdirectory(thirdparty)

# List all sources.
set(HYPERBOREAN_SUN_SOURCES
    "ssSczG-revision_submission/Buffer A.frag"
    "ssSczG-revision_submission/Buffer B.frag"
    "ssSczG-revision_submission/Buffer C.frag"
    "ssSczG-revision_submission/Common.frag"
    "ssSczG-revision_submission/Image.frag"
    "ssSczG-revision_submission/shader.ron"
)

set(MUSIC_FILE
    "${PROJECT_SOURCE_DIR}/mz - hyperborean sun mix alt 02 l2.mp3"
)

set(MIDI_FILE
    "${PROJECT_SOURCE_DIR}/hyperborean_sun_temp11.mid"
)

# Pack the shader into json format.
add_custom_target(ssSczG-revision_submission.json
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS ${HYPERBOREAN_SUN_SOURCES}
    COMMAND ${SHADERBOI} pkg "${PROJECT_SOURCE_DIR}/ssSczG-revision_submission" -y
)

# Generate the player and plain shaders.
add_custom_target(revision_submission.h
    DEPENDS ssSczG-revision_submission.json ${MUSIC_FILE} ${MIDI_FILE}
    COMMAND ${SHADERTOY_TO_C} --shader "${CMAKE_CURRENT_BINARY_DIR}/ssSczG-revision_submission.json" -v -p --multiple-shaders --midi-file ${MIDI_FILE} --music ${MUSIC_FILE}
)
add_dependencies(revision_submission.h ssSczG-revision_submission.json)

# Validate plain shaders.
set(VALIDATABLE_SHADERS
    "${CMAKE_CURRENT_BINARY_DIR}/revision_submission/buffer_a.frag"
    "${CMAKE_CURRENT_BINARY_DIR}/revision_submission/buffer_b.frag"
    "${CMAKE_CURRENT_BINARY_DIR}/revision_submission/buffer_c.frag"
    "${CMAKE_CURRENT_BINARY_DIR}/revision_submission/image.frag"
)
add_custom_target(validate-shaders
    DEPENDS revision-submission.h
)
foreach(VALIDATABLE_SHADER ${VALIDATABLE_SHADERS})
    add_custom_command(TARGET validate-shaders
        COMMAND ${GLSLANG_VALIDATOR} ${VALIDATABLE_SHADER}
    )
endforeach()
add_dependencies(validate-shaders revision_submission.h)

# Player executable
set(CMAKE_CXX_STANDARD 20 REQUIRED)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_executable(hyperborean-sun
    main.cpp
)
add_dependencies(hyperborean-sun revision_submission.h validate-shaders)
target_include_directories(hyperborean-sun PUBLIC
    "${PROJECT_SOURCE_DIR}/thirdparty/minimp3"
    "${PROJECT_SOURCE_DIR}/thirdparty/midifile/include"
)
target_link_libraries(hyperborean-sun PUBLIC
    OpenGL::GL
    CLI11
    winmm
    midifile
)
