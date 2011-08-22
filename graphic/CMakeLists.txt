
find_package(OpenGL REQUIRED)
find_package(SDL REQUIRED)
find_package(Freetype REQUIRED)
#find_package(SDL_ttf REQUIRED)
# required: FTGL
set(FTGL_LIBRARY ftgl)


#include_directories(${OPENGL_INCLUDE_DIR} ${SDL_INCLUDE_DIR} ${SDLTTF_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR} ${SDL_INCLUDE_DIR} ${FREETYPE_INCLUDE_DIRS})

set(graphic_SRCS graphic.cpp Window.cpp)

add_executable(graphic ${graphic_SRCS})
#target_link_libraries(graphic ${OPENGL_LIBRARIES} ${SDL_LIBRARY} ${SDLTTF_LIBRARY})
target_link_libraries(graphic ${OPENGL_LIBRARIES} ${SDL_LIBRARY} ${FTGL_LIBRARY})
