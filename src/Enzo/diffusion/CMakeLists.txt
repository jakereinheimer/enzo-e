# See LICENSE_CELLO file for license and copyright information

# add  to the enzo target

# create Enzo::diffusion, which is built from uncategorized source files.
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_diffusion
  EnzoMethodDiffusion.cpp EnzoMethodDiffusion.hpp
  EnzoMethodConduction.cpp
  EnzoMethodViscosity.cpp
)
add_library(Enzo::diffusion ALIAS Enzo_diffusion)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_diffusion PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_diffusion PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_diffusion PRIVATE ${Cello_TARGET_LINK_OPTIONS})
