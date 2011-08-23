#ifndef graphic_Window_h
#define graphic_Window_h

#include <iostream>
#include <cstdlib>
#include <algorithm>

#include <SDL/SDL.h>
#include <GL/gl.h>
#include <FTGL/ftgl.h>

#include "Graph.hpp"

class Graph;

namespace text {
	enum Align {
		left, right, center
	};
	enum VAlign {
		bottom, top, middle
	};
}

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
	
private:
	bool init();
	void render();
	void onLoop();
	void onEvent(SDL_Event*);
	
private:
	unsigned width, height, color;
	bool active;
	SDL_Surface *surface;
	FTTextureFont font;
	//FTPolygonFont font;
	FTSimpleLayout ftlayout;
	
	Graph *graph;
};

#endif
