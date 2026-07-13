cmake_minimum_required(VERSION 3.15)
project(user_residual CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(residual residual.cpp)

# Put the binary NEXT TO resasm.yml, not in build/: resasm.yml runs it as
# `./residual` relative to the config directory. (Windows adds the .exe itself.)
set_target_properties(residual PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY                 ${CMAKE_CURRENT_SOURCE_DIR}
  RUNTIME_OUTPUT_DIRECTORY_DEBUG           ${CMAKE_CURRENT_SOURCE_DIR}
  RUNTIME_OUTPUT_DIRECTORY_RELEASE         ${CMAKE_CURRENT_SOURCE_DIR}
  RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO  ${CMAKE_CURRENT_SOURCE_DIR}
  RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL      ${CMAKE_CURRENT_SOURCE_DIR})

# To build with OTILib sensitivities, install OTILib
# (https://github.com/mauriaristi/otilib.git, GPLv3) and add:
#   target_include_directories(residual PRIVATE ${OTILIB_INCLUDE_DIR}
#                              ${CMAKE_SOURCE_DIR}/../../partner_kit/include)
#   target_compile_definitions(residual PRIVATE RESASM_HAVE_OTILIB)
#   target_link_directories(residual PRIVATE ${OTILIB_LIBRARY_DIR})
#   target_link_libraries(residual PRIVATE oti)
