aux_source_directory(src qca_sources)

add_library (qca SHARED ${da_sources})
target_link_libraries (qca)

install (TARGETS qca LIBRARY DESTINATION lib)
