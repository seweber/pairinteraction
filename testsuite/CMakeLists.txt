# Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction library.
#
# The pairinteraction library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.

# Build tests only when testing
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL ON)

# To add test cases easier we use the following function
function(PYTHON_TEST)
  # Parse arguments
  cmake_parse_arguments(UT "" "TARGET" "SOURCE" ${ARGN} )
  # Copy test to binary directory
  configure_file(${UT_SOURCE} ${CMAKE_CURRENT_BINARY_DIR}/${UT_SOURCE})
  # Tests are executed in the root directory
  add_test(NAME "python.${UT_TARGET}" COMMAND ${PYTHON_EXECUTABLE} ${UT_SOURCE})
  set_tests_properties("python.${UT_TARGET}"
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR};${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}")
endfunction(PYTHON_TEST)

# Add Python test cases
if (WITH_PYTHON)
  add_dependencies(check _pireal _picomplex)
  find_package(PythonInterp 3 REQUIRED)
  python_test(TARGET explicitbasis SOURCE explicitbasis.py)
  python_test(TARGET greentensor SOURCE greentensor.py)
  python_test(TARGET fieldcombination SOURCE fieldcombination.py)
  python_test(TARGET feast SOURCE feast.py)
  python_test(TARGET rotation SOURCE rotation.py)
  python_test(TARGET cache SOURCE cache.py)
  python_test(TARGET perturbation SOURCE perturbation.py)
  python_test(TARGET states SOURCE state.py)
  python_test(TARGET compute SOURCE compute.py)
  python_test(TARGET quantum_defect SOURCE quantum_defect.py)
  python_test(TARGET wavefunction SOURCE wavefunction.py)
  python_test(TARGET symmetries SOURCE symmetries.py)
  python_test(TARGET integration SOURCE integration.py)
  if(NOT MSVC AND NOT (APPLE AND DEFINED ENV{TRAVIS}) AND NOT WITH_CLANG_TIDY) # timeout
    python_test(TARGET parallelization SOURCE parallelization.py)
  endif()
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/integration_test_referencedata.pickle
  ${CMAKE_CURRENT_BINARY_DIR}/integration_test_referencedata.pickle
  COPYONLY)
