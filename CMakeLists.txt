__files_group(${CMAKE_CURRENT_SOURCE_DIR}/cxutil SRC)

__cc_find(Clipper)
__cc_find(RapidJson)

__assert_target(clipper)
#__assert_parameter(RAPIDJSON_INCLUDE_DIRS)

set(LIBS clipper trimesh2 stringutil rapidjson)
set(DEFS _CRT_SECURE_NO_WARNINGS NOMINMAX)
set(INCS )
__add_real_target(cxutil lib SOURCE ${SRC}
							 LIB ${LIBS}
							 DEF ${DEFS}
							 INC ${INCS}
							 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							 OPENMP
							 )
							 
__set_folder_targets(slice TARGET cxutil)



