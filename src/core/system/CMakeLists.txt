add_library(SF_formula STATIC SF_formula.cpp SF_formula.h SF_provenance.h)
target_link_libraries(SF_formula PUBLIC SF_headers)
set_lib_output(SF_formula core)
