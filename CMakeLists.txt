cmake_minimum_required(VERSION 3.31.1)

# 修改项目名称，避免特殊字符造成解析错误
project(matplot_test)

include_directories(/usr/local/include/eigen-3.4.0)


set(CMAKE_CXX_STANDARD 17)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

find_package(Eigen3 3.4 REQUIRED)
if(TARGET Eigen3::Eigen)
message(STATUS "*** Eigen3 found ***")
endif()

# 查找 Matplot++ 库
find_package(Matplot++ REQUIRED)

# 如果找到 Matplot++，打印信息
if (TARGET Matplot++::matplot)
    message(STATUS "*** Matplot++ found ***")
else()
    message(FATAL_ERROR "Matplot++ not found!")
endif()

# 创建可执行文件
add_executable(matplottest main.cc)

# 链接 Matplot++ 库
target_link_libraries(matplottest
    PUBLIC
    Matplot++::matplot
)