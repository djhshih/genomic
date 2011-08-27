find_package(OpenGL REQUIRED)
find_package(SDL REQUIRED)
find_package(Freetype REQUIRED)

find_package(FTGL REQUIRED)


include_directories(${OPENGL_INCLUDE_DIR} ${SDL_INCLUDE_DIR} ${FREETYPE_INCLUDE_DIRS} ${FTGL_INCLUDE_DIR})

set(graphic_SRCS graphic.cpp Window.cpp Graph.cpp)

add_executable(graphic ${graphic_SRCS})
target_link_libraries(graphic ${OPENGL_LIBRARIES} ${SDL_LIBRARY} ${FTGL_LIBRARY})
