#ifndef command_graphic_h
#define command_graphic_h

#include "../graphic/Window.hpp"
#include "../graphic/Graph.hpp"

extern Window window;

// tolua_begin

//typedef Graph::x_type x_type;
//typedef Graph::y_type y_type;
typedef long x_type;
typedef float y_type;

void jumpTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest);

void moveTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest);

void easeTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest);

// tolua_end


#endif
