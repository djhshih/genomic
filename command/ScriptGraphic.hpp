#ifndef command_ScriptGraphic_h
#define command_ScriptGraphic_h

#include "../graphic/Window.hpp"
#include "Script.hpp"
#include "graphic.hpp"

int tolua_graphic_open (lua_State* tolua_S);

class ScriptGraphic : public Script
{
public:
	
	ScriptGraphic(Window& _window)
	: window(_window) {
		tolua_graphic_open(L);	
	}
	
private:
	
	Window& window;
	
};

#endif
