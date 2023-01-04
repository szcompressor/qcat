file(GLOB cpu_source_files "${CMAKE_SOURCE_DIR}/examples/*.c")

foreach (SRCFILE IN LISTS cpu_source_files)
    get_filename_component(EXE ${SRCFILE} NAME_WLE)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.c)
    target_link_libraries(${EXE} PRIVATE qcat)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${EXE}$<TARGET_FILE_SUFFIX:${EXE}> TYPE BIN)

endforeach ()