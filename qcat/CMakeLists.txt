aux_source_directory(src qcat_sources)

add_library (qcat SHARED ${da_sources})
target_link_libraries (qcat)

install (TARGETS qcat LIBRARY DESTINATION lib)
