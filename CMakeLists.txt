__files_group(${CMAKE_CURRENT_SOURCE_DIR}/cxutil SRC)

__cc_find(clipper)
__cc_find(rapidjson)

__assert_target(clipper)
__assert_target(rapidjson)

set(LIBS clipper trimesh2 rapidjson)
set(DEFS _CRT_SECURE_NO_WARNINGS NOMINMAX)
set(INCS )
__add_real_target(cxutil lib SOURCE ${SRC}
							 ILIB ${LIBS}
							 DEF ${DEFS}
							 INC ${INCS}
							 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							 OPENMP
							 )
							 



