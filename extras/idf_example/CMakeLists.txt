# ESP-IDF example: builds the BasicRender benchmark as an IDF app, using the FreqMul component two levels up.
#   cd extras/idf_example && idf.py set-target esp32s3 && idf.py build flash monitor
cmake_minimum_required(VERSION 3.16)
list(APPEND EXTRA_COMPONENT_DIRS "${CMAKE_CURRENT_LIST_DIR}/../..")
include($ENV{IDF_PATH}/tools/cmake/project.cmake)
project(freqmul_idf_example)
