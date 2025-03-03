include(TestSetup)

set(_core_libs IIsrc IIparallel IItools IIxdr)

set(TESTS_NAMES "")

# AtomSelection tests
set(AtomSelectionTest "")

# Finding the git package
find_package(Git)
if(Git_FOUND)
  message("Git found: ${GIT_EXECUTABLE}")
endif()

add_unit_test("TestRotation"
              "TestRotation"
              SOURCES main.cpp
              LIBS ${_core_libs}
              ARGS ${CMAKE_CURRENT_LIST_DIR}/TestRotation/input.dat
)

# add_unit_test("TestAtomSelection/TestIndexFileSelection"
#               "TestIndexSelection"
#               SOURCES main.cpp
#               LIBS ${_core_libs}
#               ARGS ${CMAKE_CURRENT_LIST_DIR}/TestAtomSelection/TestIndexFileSelection/input.dat
# )

# Now add the test which tests the main program (only test if they run) 
set(OMP_THREADS 1 4 8)
set(OP_TEST "")

set(PROGRAM_PATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/mdanalysis")
set(TEST_RUN_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/run_test.sh")

# Add the test and run with the bash script run_test.sh
foreach(TEST_NAME IN LISTS MAIN_TESTS)
    foreach(THREADS IN LISTS OMP_THREADS)
        set(TEST_SCRIPT "${TEST_NAME}/input.dat")
        set(TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/${TEST_NAME}")
        set(TEST_REF_DIR "${CMAKE_CURRENT_LIST_DIR}/SampleTrajectory")
        set(TEST_CMD ${BASH} ${TEST_RUN_SCRIPT} -p ${PROGRAM_PATH} -d ${TEST_DIR} 
        -o ${CMAKE_CURRENT_BINARY_DIR}/op.out -r ref.out
             -i ${CMAKE_CURRENT_LIST_DIR}/${TEST_SCRIPT} -n ${THREADS} -f ${TEST_REF_DIR})

        # set working directory to be the directory processed currently
        add_test(NAME "${TEST_NAME}::OMP_${THREADS}"
                COMMAND ${TEST_CMD}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
        ) 
    endforeach()
endforeach()

add_custom_target(build_test COMMAND ${CMKAE_CTEST_COMMAND}
                              DEPENDS ${TESTS_NAMES})