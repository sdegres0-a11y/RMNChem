cmake_minimum_required(VERSION 3.15)
project(RMNChem C)

include(GenerateG3A)
include(Fxconv)
find_package(Gint 2.9 REQUIRED)

set(SOURCES
  src/main.c
)

add_executable(rmnchem ${SOURCES})

# -Os plutot que -O2 : l'add-in est limite par sa taille, pas par le CPU.
# Le decoupage en sections + --gc-sections elimine tout ce qui n'est pas
# reellement appele, ce qui reduit sensiblement le binaire final.
target_compile_options(rmnchem PRIVATE
  -Wall -Wextra -Os
  -ffunction-sections -fdata-sections
)
target_link_options(rmnchem PRIVATE -Wl,--gc-sections)

target_link_libraries(rmnchem Gint::Gint)

if("${FXSDK_PLATFORM_LONG}" STREQUAL fxCG50)
  generate_g3a(TARGET rmnchem OUTPUT "RMNChem.g3a"
    NAME "RMNChem"
    ICONS assets-cg/icon-uns.png assets-cg/icon-sel.png)
endif()
