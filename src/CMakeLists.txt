__files_group(${CMAKE_CURRENT_SOURCE_DIR}/cxutil SRC)

__assert_target(clipper)
__assert_target(spdlog)

set(LIBS clipper spdlog)
set(DEFS _CRT_SECURE_NO_WARNINGS NOMINMAX)

__add_real_target(cxutil lib SOURCE ${SRC} LIB ${LIBS} DEF ${DEFS})
set_property(TARGET cxutil PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})



