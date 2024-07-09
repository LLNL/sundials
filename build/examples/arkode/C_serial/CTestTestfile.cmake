# CMake generated Testfile for 
# Source directory: /home/maggul/STS/sundials/examples/arkode/C_serial
# Build directory: /home/maggul/STS/sundials/build/examples/arkode/C_serial
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ark_analytic "/home/maggul/STS/sundials/build/examples/arkode/C_serial/ark_analytic")
set_tests_properties(ark_analytic PROPERTIES  _BACKTRACE_TRIPLES "/home/maggul/STS/sundials/cmake/macros/SundialsAddTest.cmake;204;add_test;/home/maggul/STS/sundials/examples/arkode/C_serial/CMakeLists.txt;152;sundials_add_test;/home/maggul/STS/sundials/examples/arkode/C_serial/CMakeLists.txt;0;")
add_test(lsrk_analytic "/home/maggul/STS/sundials/build/examples/arkode/C_serial/lsrk_analytic")
set_tests_properties(lsrk_analytic PROPERTIES  _BACKTRACE_TRIPLES "/home/maggul/STS/sundials/cmake/macros/SundialsAddTest.cmake;204;add_test;/home/maggul/STS/sundials/examples/arkode/C_serial/CMakeLists.txt;152;sundials_add_test;/home/maggul/STS/sundials/examples/arkode/C_serial/CMakeLists.txt;0;")
add_test(lsrk_analytic_VarJac "/home/maggul/STS/sundials/build/examples/arkode/C_serial/lsrk_analytic_VarJac")
set_tests_properties(lsrk_analytic_VarJac PROPERTIES  _BACKTRACE_TRIPLES "/home/maggul/STS/sundials/cmake/macros/SundialsAddTest.cmake;204;add_test;/home/maggul/STS/sundials/examples/arkode/C_serial/CMakeLists.txt;152;sundials_add_test;/home/maggul/STS/sundials/examples/arkode/C_serial/CMakeLists.txt;0;")
add_test(erk_analytic "/home/maggul/STS/sundials/build/examples/arkode/C_serial/erk_analytic")
set_tests_properties(erk_analytic PROPERTIES  _BACKTRACE_TRIPLES "/home/maggul/STS/sundials/cmake/macros/SundialsAddTest.cmake;204;add_test;/home/maggul/STS/sundials/examples/arkode/C_serial/CMakeLists.txt;152;sundials_add_test;/home/maggul/STS/sundials/examples/arkode/C_serial/CMakeLists.txt;0;")
