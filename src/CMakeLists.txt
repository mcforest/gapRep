# add_executable (hello hallo.cpp)


add_executable (precalc precalculation.cpp)
target_link_libraries    (precalc -ldl)
target_link_libraries    (precalc glog ${GLOG_LIBRARY})
target_link_libraries    (precalc gflags ${GFLAGS_LIBRARY})
target_link_libraries    (precalc ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries    (precalc sdsl)


add_custom_command (
    TARGET precalc
    COMMENT "bestimme grenzwerte..."
    COMMAND precalc
    DEPENDS precalc
)


add_executable (gapRep2 main.cpp)
target_link_libraries    (gapRep2 -ldl)
target_link_libraries    (gapRep2 glog ${GLOG_LIBRARY})
target_link_libraries    (gapRep2 gflags ${GFLAGS_LIBRARY})
target_link_libraries    (gapRep2 ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries    (gapRep2 sdsl)

add_dependencies (gapRep2 precalc)


