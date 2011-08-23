#ifndef graphic_Graph_h
#define graphic_Graph_h

#include <cstdio>
#include <vector>
#include <stdexcept>

#include <GL/gl.h>

#include "Window.hpp"

class Window;

class Graph
{
public:
	typedef long xtype;
	typedef float ytype;
	
	float width, height;
	
public:
	
	Graph(Window &parentWindow)
	: parent(parentWindow), nlists(1),
	  xmin(0), xmax(10000), xstep(1000),
	  ymin(-2), ymax(4), ystep(1)
	{
		init();
	}
	
	void plot(const std::vector<xtype>& x, const std::vector<ytype>& y) {
		if (x.size() == 0) return;
		if (x.size() != y.size()) {
			throw std::runtime_error("x and y cannot differ in size");
		}
		glPointSize(2);
		float rangexf = float(xmax-xmin);
		float rangeyf = float(ymax-ymin);
		glBegin(GL_POINTS);
		for (std::size_t i = 0; i < x.size(); ++i) {
			glVertex2f(
				(x[i]-xmin)/rangexf,
				(y[i]-ymin)/rangeyf
			);
		}
		glEnd();
		
		glColor3f(0, 1, 0);
		glPointSize(1);
		glBegin(GL_LINES);
		for (std::size_t i = 0; i < x.size(); ++i) {
			glVertex2f(
				(x[i]-xmin)/rangexf,
				(y[i]-ymin)/rangeyf
			);
		}
		glEnd();
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
	/*
	ytype rescale(ytype y) {
		return (y-ymin) * 1.0f/ymax + ymin/ymax;
	}
	*/
	
private:
	
	Window &parent;
	GLuint lists, nlists;
	xtype xmin, xmax, xstep;
	ytype ymin, ymax, ystep;
	
};

#endif
