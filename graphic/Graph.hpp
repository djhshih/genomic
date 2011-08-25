#ifndef graphic_Graph_h
#define graphic_Graph_h

#include <cstdio>
#include <vector>
#include <stdexcept>
#include <queue>

#include <SDL/SDL.h>
#include <GL/gl.h>

#include "Window.hpp"
#include "Tween.hpp"

class Window;


class Graph
{
public:
	typedef long x_type;
	typedef float y_type;
	
	struct PropertiesDraw {
		bool xaxis;
		bool yaxis;
		bool data_points;
		bool data_lines; 
		bool background;
		bool reference_line;
		PropertiesDraw()
		: xaxis(true), yaxis(true), background(true), data_points(true), data_lines(true), reference_line(true) {}
	};
	
	struct Area {
		x_type xstart, xend;
		y_type ystart, yend;
		Area(x_type xb, x_type xe, y_type yb, y_type ye)
		: xstart(xb), xend(xe), ystart(yb), yend(ye) {}
	};
	
	typedef bool (Graph::*task)();
	typedef std::queue<task> taskqueue;
	typedef std::queue<Area> destqueue;
	
	float width, height;
	
public:
	
	Graph(Window &parentWindow)
	: parent(parentWindow), nlists(1),
	  xmin(0), xmax(10000),
	  ymin(-2), ymax(4),
	  xstart(0), xend(10000), xstep(1), xintervals(10),
	  ystart(-2), yend(4), ystep(0.5), yintervals(6),
	  xtickmajor(500), ytickmajor(1),
	  tickSize(0.01),
	  labelOffset(0.03), titleOffset(0.12), axisOffset(0.05),
	  bgOverflow(0.02),
	  xease(20), yease(20),
	  complete(true),
	  delay(10),
	  duration(100),
	  width(1), height(.5)
	{
		init();
	}
	
	void scrollX() {
		xstart += xstep;
		xend += xstep;
		if (xstart < xmin) xstart = xmin;
		if (xend > xmax) xend = xmax;
		SDL_Delay(delay);
	}
	
	void scrollY() {
		ystart += ystep;
		yend += ystep;
		if (ystart < ymin) ystart = ymin;
		if (yend > ymax) yend = ymax;
		SDL_Delay(delay);
	}
	
	void addTask(task t, x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
		dests.push(Area(xstartDest, xendDest, ystartDest, yendDest));
		tasks.push(t);
	}
	
	void jumpTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
		addTask(&Graph::jumpToTask, xstartDest, xendDest, ystartDest, yendDest);
	}
	
	bool jumpToTask() {
		xstart += xstart_change;
		xend += xend_change;
		ystart += ystart_change;
		yend += yend_change;
		if (xstart < xmin) xstart = xmin;
		if (ystart < ymin) ystart = ymin;
		if (xend > xmax) xend = xmax;
		if (yend > ymax) yend = ymax;
		return true;
	}
	
	void moveTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
		addTask(&Graph::moveToTask, xstartDest, xendDest, ystartDest, yendDest);
	}
	
	bool moveToTask() {
		xstart = tween::linear(time, xstart_src, xstart_change, duration);
		xend = tween::linear(time, xend_src, xend_change, duration);
		ystart = tween::linear(time, ystart_src, ystart_change, duration);
		yend = tween::linear(time, yend_src, yend_change, duration);
		if (++time > duration) {
			return true;
		}
		return false;
	}
	
	void easeTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
		addTask(&Graph::easeToTask, xstartDest, xendDest, ystartDest, yendDest);
	}
	
	bool easeToTask() {
		xstart = tween::easeCubic(time, xstart_src, xstart_change, duration);
		xend = tween::easeCubic(time, xend_src, xend_change, duration);
		ystart = tween::easeCubic(time, ystart_src, ystart_change, duration);
		yend = tween::easeCubic(time, yend_src, yend_change, duration);
		if (++time > duration) {
			return true;
		}
		return false;
	}
	
	void update() {
		// member function pointer must be calle using (obj.*p)(...)
		if (!tasks.empty()) {
			pre();
			if ( ( this->*(tasks.front()) )() ) {
				// function returns success: remove it from list
				tasks.pop();
				std::cout << "done\n";
				post();
			}
		}
	}
	
	/*
	void update() {
		taskqueue::iterator it = tasks.start();
		taskqueue::const_iterator end = tasks.end();
		while (it != end) {
			// call function in list
			// dereference iterator to get member function pointer, which must be calle using (obj.*p)(...)
			if ( ( this->*(*it) )() ) {
				// function returns success: remove it from list
				it = tasks.erase(it);
				// post-processing
				post();
			} else {
				++it;
			}
		}
	}
	*/
	
	void pre() {
		if (complete) {
			Area dest = dests.front();
			dests.pop();
			xstart_change = dest.xstart - xstart;
			xend_change = dest.xend - xend;
			ystart_change = dest.ystart - ystart;
			yend_change = dest.yend - yend;
			
			xstart_src = xstart;
			xend_src = xend;
			ystart_src = ystart;
			yend_src = yend;
			
			time = 0;
			complete = false;
		}
	}
	
	void post() {
		xtickmajor = (xend - xstart) / xintervals;
		ytickmajor = (yend - ystart) / yintervals;
		complete = true;
	}
	
	/*
	void easeto_x(x_type start, x_type end) {
		if (start < xmin) start = xmin;
		if (end > xmax) end = xmax;
		
		xstart += (start - xstart) / xease;
		xend += (end - xend) / xease;
		xtickmajor = (xend - xstart) / xintervals;
		
		if (abs(xstart - start) < xease) {
			xstart = start;
			xend = end;
			draws.xaxis = true;
		} else {
			draws.xaxis = false;
		}
		
		SDL_Delay(delay);
	}
	*/
	
	void plot(const std::vector<x_type>& x, const std::vector<y_type>& y);
	
	void drawTicksX();
	void drawTicksY();
	
	GLuint begin() {
		return lists;
	}
	
	GLuint end() {
		return lists + nlists;
	}
	
private:
	
	void init();
	
	void jumpToX(x_type start, x_type end) {
		if (start >= xmin) xstart = start;
		if (end <= xmax) xend = end;
		post();
	}
	
	void jumpToY(y_type start, y_type end) {
		if (start >= ymin) ystart = start;
		if (end <= ymax) yend = end;
		post();
	}
	
private:
	
	Window &parent;
	GLuint lists, nlists;
	x_type xstart, xend, xstep, xtickmajor, xease;
	x_type xstart_src, xstart_change, xend_src, xend_change;
	const x_type xmin, xmax, xintervals;
	y_type ystart, yend, ystep, ytickmajor, yease;
	y_type ystart_src, ystart_change, yend_src, yend_change;
	const y_type ymin, ymax, yintervals;
	const float tickSize, labelOffset, titleOffset, axisOffset, bgOverflow;
	const unsigned delay;
	float time;
	const float duration;
	
	taskqueue tasks;
	destqueue dests;
	
	bool complete;
	
	PropertiesDraw draws;
	
};

#endif
