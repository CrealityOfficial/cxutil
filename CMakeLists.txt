__files_group(${CMAKE_CURRENT_SOURCE_DIR}/cxutil SRC)

__assert_target(clipper)
__assert_target(spdlog)
__assert_parameter(RAPIDJSON_INCLUDE_DIRS)

set(LIBS clipper spdlog trimesh2)
if(TARGET OpenMP::OpenMP_CXX)
	list(APPEND LIBS OpenMP::OpenMP_CXX)
endif()

set(DEFS _CRT_SECURE_NO_WARNINGS NOMINMAX)
set(INCS ${RAPIDJSON_INCLUDE_DIRS})
__add_real_target(cxutil lib SOURCE ${SRC}
							 LIB ${LIBS}
							 DEF ${DEFS}
							 INC ${INCS}
							 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							 )



