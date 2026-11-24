#--------------------------------------------------------------------------
# Directory Setup
#--------------------------------------------------------------------------
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

#--------------------------------------------------------------------------
# Sphinx Configuration
#--------------------------------------------------------------------------
find_package(Sphinx REQUIRED)

# Main Sphinx documentation build
add_custom_target(sphinx
  COMMAND ${SPHINX_EXECUTABLE} -M html
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  -c ${SPHINX_SOURCE} -W --keep-going -n
  COMMENT "Generating HTML documentation with Sphinx"
)
#--------------------------------------------------------------------------
# Main Documentation Target
#--------------------------------------------------------------------------
add_custom_target(docs ALL)
add_dependencies(docs sphinx)
