#ifndef graphic_Graph_h
#define graphic_Graph_h

#include <GL/gl.h>
#include "Window.hpp"

class Window;

class Graph
{
public:
	float width, height;
	
public:
	
	Graph(Window &parentWindow)
	: parent(parentWindow), nlists(1),
	  xmin(0), xmax(10000), xstep(1000),
	  ymin(-2), ymax(4), ystep(1)
	{
		init();
	}
	
	GLuint begin() {
		return lists;
	}
	GLuint end() {
		return lists + nlists;
	}
	
private:
	
	void init();
	
	// rescale y s.t. distance(0, ymax) == 1
	float rescale(float y) {
		return (y-ymin) * 1.0f/ymax + ymin/ymax;
	}
	
private:
	
	Window &parent;
	GLuint lists, nlists;
	long xmin, xmax, xstep;
	float ymin, ymax, ystep;
	
};

#endif
