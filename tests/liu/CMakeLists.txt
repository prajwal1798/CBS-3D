cmake_minimum_required(VERSION 3.16)
project(CBS3D_Liu_Element_Verification LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

get_filename_component(CBS3D_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)

add_executable(liu_nithiarasu_element_verification
    "${CBS3D_ROOT}/tests/LiuNithiarasuElementVerification.cpp"
    "${CBS3D_ROOT}/src/assembly/MomentumAssemblyLiu.cpp"
    "${CBS3D_ROOT}/src/assembly/SpalartAllmarasAssemblyLiu.cpp"
    "${CBS3D_ROOT}/src/turbulence/SpalartAllmaras.cpp"
)

target_include_directories(liu_nithiarasu_element_verification PRIVATE
    "${CBS3D_ROOT}/include"
)

if (MSVC)
    target_compile_options(liu_nithiarasu_element_verification PRIVATE
        /W4
        /WX
        /permissive-
    )
else()
    target_compile_options(liu_nithiarasu_element_verification PRIVATE
        -O2
        -Wall
        -Wextra
        -Werror
        -pedantic
    )
endif()

enable_testing()
add_test(
    NAME liu_nithiarasu_element_formulation
    COMMAND liu_nithiarasu_element_verification
)
