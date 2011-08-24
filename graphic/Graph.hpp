#ifndef graphic_Graph_h
#define graphic_Graph_h

#include <cstdio>
#include <vector>
#include <stdexcept>

#include <SDL/SDL.h>
#include <GL/gl.h>

#include "Window.hpp"

class Window;

class Graph
{
public:
	typedef long x_type;
	typedef float y_type;
	
	float width, height;
	
	struct draw {
		enum drawflag {
			xaxis = 0x01,
			yaxis = 0x02,
			background = 0x04
		};
	};
	
public:
	
	Graph(Window &parentWindow)
	: parent(parentWindow), nlists(1),
	  xmin(0), xmax(10000),
	  ymin(-2), ymax(4),
	  xbegin(0), xend(10000), xstep(1), xintervals(10),
	  ybegin(-2), yend(4), ystep(0.5), yintervals(6),
	  xtickmajor(500), ytickmajor(1),
	  tickSize(0.01),
	  labelOffset(0.03), titleOffset(0.12), axisOffset(0.05),
	  bgOverflow(0.02),
	  xease(20), yease(20),
	  delay(10),
	  skipflags(0),
	  width(1), height(.5)
	{
		init();
	}
	
	void scroll_x() {
		xbegin += xstep;
		xend += xstep;
		SDL_Delay(delay);
	}
	
	void scroll_y() {
		ybegin += ystep;
		yend += ystep;
		SDL_Delay(delay);
	}
	
	void goto_x(x_type begin, x_type end) {
		if (begin >= xmin) xbegin = begin;
		if (end <= xmax) xend = end;
		xtickmajor = (xend - xbegin) / xintervals;
	}
	
	void goto_y(y_type begin, y_type end) {
		if (begin >= ymin) ybegin = begin;
		if (end <= ymax) yend = end;
		ytickmajor = (yend - ybegin) / yintervals;
	}
	
	void easeto_x(x_type begin, x_type end) {
		if (begin < xmin) begin = xmin;
		if (end > xmax) end = xmax;
		
		xbegin += (begin - xbegin) / xease;
		xend += (end - xend) / xease;
		xtickmajor = (xend - xbegin) / xintervals;
		
		if (abs(xbegin - begin) < xease) {
			skipflags &= ~(draw::xaxis);
		} else {
			skipflags |= draw::xaxis;
		}
		
		SDL_Delay(delay);
	}
	
	void plot(const std::vector<x_type>& x, const std::vector<y_type>& y);
	
	GLuint begin() {
		return lists;
	}
	
	GLuint end() {
		return lists + nlists;
	}
	
private:
	
	void init();
	
private:
	
	Window &parent;
	GLuint lists, nlists;
	x_type xbegin, xend, xstep, xtickmajor, xease;
	const x_type xmin, xmax, xintervals;
	y_type ybegin, yend, ystep, ytickmajor, yease;
	const y_type ymin, ymax, yintervals;
	const float tickSize, labelOffset, titleOffset, axisOffset, bgOverflow;
	const unsigned delay;
	
	int skipflags;
	
};

#endif
