__files_group(${CMAKE_CURRENT_SOURCE_DIR}/cxutil SRC)

__cc_find(Clipper)
__cc_find(TriMesh2)
__cc_find(RapidJson)

__assert_target(clipper)
__assert_target(trimesh2)
__assert_target(rapidjson)

set(ILIBS clipper trimesh2)
if(TARGET rapidjson)
	list(APPEND ILIBS rapidjson)
endif()

set(DEFS _CRT_SECURE_NO_WARNINGS NOMINMAX)
set(INCS ${RAPIDJSON_INCLUDE_DIRS})
__add_real_target(cxutil lib SOURCE ${SRC}
							 ILIB ${ILIBS}
							 DEF ${DEFS}
							 INC ${INCS}
							 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							 OPENMP
							 )
							 
__set_folder_targets(slice TARGET cxutil)



