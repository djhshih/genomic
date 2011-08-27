find_package(Lua51 REQUIRED)

include_directories(${LUA_INCLUDE_DIRS})

set(command_SRCS command.cpp Script.cpp)

add_executable(command ${command_SRCS})
target_link_libraries(command ${LUA_LIBRARIES})
