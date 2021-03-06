#ifndef graphic_Window_h
#define graphic_Window_h

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <stdexcept>
#include <vector>

#include <SDL/SDL.h>
#include <SDL/SDL_thread.h>
#include <GL/gl.h>
#include <FTGL/ftgl.h>

#include "Graph.hpp"

namespace text {
	enum Align {
		left, right, center
	};
	enum VAlign {
		bottom, top, middle
	};
}

//TODO fix circular dependency of Graph and Window
//TODO   which appears to cause problems with delete graph

class Graph;

class Window
{
public:
	
	typedef int (*thread_function)(void *);
	
	Window()
	: surface(NULL), active(true),
	  width(1024), height(748), color(32),
	  font("/usr/share/fonts/TTF/Vera.ttf"),
	  graph(NULL)
	{
		if (!init()) {
			throw std::runtime_error("Failed to initialize window.");
		}
	}
	
	~Window();
	
	int exec();
	
	void addThread(thread_function function) {
		if (function != NULL) {
			SDL_Thread *thread = SDL_CreateThread(function, NULL);
			if (thread == NULL) {
				throw std::runtime_error("Failed to create thread.");
			}
			threads.push_back(thread);
		}
	}
	
	void renderText(const char *text, float x=0, float y=0, float angle=0, unsigned int size=10, text::Align align=text::center, text::VAlign valign=text::middle) {
		glPushMatrix();
		
		// perform absolute translation before rotation
		glTranslatef(x, y, 0);
		
		// scale the font to the coordinate system
		glScalef(1.0f/width, 1.0f/height, 1);
		font.FaceSize(size);
		
		// perform rotation
		glRotatef(angle, 0, 0, 1);
		
		// perform relative translation after rotation
		FTBBox box;
		if (align != text::left || valign != text::bottom) {
			box = font.BBox(text);
		}
		switch (align) {
			case text::center:
				glTranslatef(-(box.Upper().Xf() - box.Lower().Xf())/2, 0, 0);
				break;
			case text::right:
				glTranslatef(-(box.Upper().Xf() - box.Lower().Xf()), 0, 0);
				break;
			default:
				break;
		}
		switch (valign) {
			case text::middle:
				glTranslatef(0, -(box.Upper().Yf() - box.Lower().Yf())/2, 0);
				break;
			case text::top:
				glTranslatef(0, -(box.Upper().Yf() - box.Lower().Yf()), 0);
				break;
			default:
				break;
		}
		
		font.Render(text);
		
		glPopMatrix();
	}
	
	void renderTextbox(const char *text, unsigned int size=10) {
		glPushMatrix();
		glScalef(0.001, 0.001, 1);
		font.FaceSize(size);
		ftlayout.Render(text);
		glPopMatrix();
	}
	
public:
	
	Graph *graph;
	
private:
	
	bool init();
	void render();
	virtual void onLoop();
	void onEvent(SDL_Event*);
	
private:
	
	unsigned width, height, color;
	bool active;
	SDL_Surface *surface;
	FTTextureFont font;
	FTSimpleLayout ftlayout;
	
	std::vector<SDL_Thread*> threads;
	
};

#endif
