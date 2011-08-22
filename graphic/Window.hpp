#ifndef graphic_Window_h
#define graphic_Window_h

#include <iostream>

#include <SDL/SDL.h>
#include <GL/gl.h>
#include <FTGL/ftgl.h>

#include "Graph.hpp"

class Graph;

class Window
{
public:
	
	Window()
	: surface(NULL), active(true),
	  width(640), height(480), color(32),
	  font("/usr/share/fonts/TTF/Vera.ttf"),
	  graph(NULL)
	{}
	
	~Window() {
		SDL_FreeSurface(surface);
		SDL_Quit();
		delete graph;
	}
	
	int exec();
	void renderText(const char *text, unsigned int size=10) {
		font.FaceSize(size);
		font.Render(text);
	}
	
private:
	bool init();
	void render();
	void onLoop();
	void onEvent(SDL_Event*);
	
private:
	unsigned width, height, color;
	bool active;
	SDL_Surface *surface;
	FTGLPixmapFont font;
	
	Graph *graph;
};

#endif
