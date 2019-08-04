########################################################################
IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    MESSAGE(FATAL_ERROR "Prevented in-tree built. This is bad practice.")
ENDIF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

########################################################################
# Project setup
########################################################################



include_directories(${CMAKE_CURRENT_SOURCE_DIR})

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} CURVE_FIT_SRCS)



ADD_LIBRARY(curve_fit SHARED ${CURVE_FIT_SRCS})
TARGET_LINK_LIBRARIES(curve_fit ${LIBS})
ADD_LIBRARY(curve_fit_static STATIC ${CURVE_FIT_SRCS})
set_target_properties(curve_fit_static PROPERTIES OUTPUT_NAME "curve_fit")
set_target_properties(curve_fit PROPERTIES CLEAN_DIRECT_OUTPUT 1)
set_target_properties(curve_fit_static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
set_target_properties(curve_fit PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set_target_properties(curve_fit_static PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

#IF(UNIX)
#    SET_PROPERTY(TARGET curve_fit PROPERTY VERSION ${ZYNQ_HARDWARE_VERSION})
#    SET_PROPERTY(TARGET curve_fit PROPERTY SOVERSION ${ZYNQ_HARDWARE_VERSION_MAJOR})
#ENDIF(UNIX)
#INSTALL(TARGETS curve_fit DESTINATION lib)
#INSTALL(DIRECTORY . DESTINATION include FILES_MATCHING PATTERN "*.h")
#INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/curve_fit.h DESTINATION include)
#INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/stream_func.h DESTINATION include)
#INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/record_media.h DESTINATION include)
#INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/discovery/broadcast_client.h DESTINATION include)