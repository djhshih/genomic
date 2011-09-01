find_package(Lua51 REQUIRED)
find_package(Freetype REQUIRED)
find_package(Boost REQUIRED)

include_directories(${LUA_INCLUDE_DIRS} ${FREETYPE_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS})

set(command_SRCS Script.cpp)

add_executable(command_demo command.cpp graphic.cpp graphic_lua.cpp ${command_SRCS})
target_link_libraries(command_demo ${LUA_LIBRARIES} graphic boost_thread tolua++)
