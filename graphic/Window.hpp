#ifndef graphic_Window_h
#define graphic_Window_h

#include <iostream>

#include <SDL/SDL.h>
//#include <SDL/SDL_ttf.h>
#include <GL/gl.h>
//#include <GL/glu.h>
#include <FTGL/ftgl.h>

class Window
{
public:
	
	Window() : surface(NULL), active(true), width(640), height(480), color(32), font("/usr/share/fonts/TTF/Vera.ttf") {}
	
	~Window() {
		SDL_FreeSurface(surface);
		SDL_Quit();
		//TTF_Quit();
	}
	
	int exec();
	
private:
	bool init();
	void render();
	void onLoop();
	void onEvent(SDL_Event*);
	
private:
	unsigned width, height, color;
	bool active;
	SDL_Surface *surface;
	//TTF_Font *font;
	FTGLPixmapFont font;
};

#endif
