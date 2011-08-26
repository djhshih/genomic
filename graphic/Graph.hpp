#ifndef graphic_Graph_h
#define graphic_Graph_h

#include <cstdio>
#include <vector>
#include <list>
#include <queue>
#include <stdexcept>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <SDL/SDL.h>
#include <GL/gl.h>

#include "Window.hpp"
#include "Area.hpp"
#include "Tween.hpp"


//TODO improve dynamic behaviour of axis

class Window;

class Graph
{
public:
	
	typedef long x_type;
	typedef float y_type;
	typedef Area<x_type, y_type> GraphArea;
	
	typedef bool (Graph::*task)();
	typedef boost::function<void(void)> plotter;
	
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
	
	float width, height;
	
	PropertiesDraw draws;
	
public:
	
	Graph(Window &parentWindow)
	: parent(parentWindow), nlists(1),
	  xmin(0), xmax(10000),
	  ymin(-2), ymax(4),
	  xstart(0), xend(10000), xstep(1), xintervals(10),
	  ystart(-2), yend(4), ystep(0.5), yintervals(6),
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
	
	void jumpTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
		addTask(&Graph::jumpToTask, xstartDest, xendDest, ystartDest, yendDest);
	}
	
	void moveTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
		addTask(&Graph::moveToTask, xstartDest, xendDest, ystartDest, yendDest);
	}
	
	void easeTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
		addTask(&Graph::easeToTask, xstartDest, xendDest, ystartDest, yendDest);
	}
	
	void update() {
		// member function pointer must be called using (obj.*p)(...)
		if (!tasks.empty()) {
			pre();
			if ( ( this->*(tasks.front()) )() ) {
				// function returns success: remove it from list
				tasks.pop();
				post();
			}
		}
	}
	
	void pre() {
		if (complete) {
			GraphArea dest = dests.front();
			dests.pop();
			
			src = GraphArea(xstart, xend, ystart, yend);
			change = dest - src;
			
			time = 0;
			complete = false;
		}
	}
	
	void post() {
		xtickmajor = (xend - xstart) / xintervals;
		ytickmajor = (yend - ystart) / yintervals;
		complete = true;
	}
	
	void plot();
	
	void add(const std::vector<x_type>& x, const std::vector<y_type>& y);
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
	
	void addTask(task t, x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
		dests.push(GraphArea(xstartDest, xendDest, ystartDest, yendDest));
		tasks.push(t);
	}
	
	bool jumpToTask() {
		xstart += change.xstart;
		xend += change.xend;
		ystart += change.ystart;
		yend += change.yend;
		if (xstart < xmin) xstart = xmin;
		if (ystart < ymin) ystart = ymin;
		if (xend > xmax) xend = xmax;
		if (yend > ymax) yend = ymax;
		return true;
	}
	
	bool moveToTask() {
		xstart = tween::linear(time, src.xstart, change.xstart, duration);
		xend = tween::linear(time, src.xend, change.xend, duration);
		ystart = tween::linear(time, src.ystart, change.ystart, duration);
		yend = tween::linear(time, src.yend, change.yend, duration);
		if (++time > duration) {
			return true;
		}
		return false;
	}
	
	bool easeToTask() {
		xstart = tween::easeCubic(time, src.xstart, change.xstart, duration);
		xend = tween::easeCubic(time, src.xend, change.xend, duration);
		ystart = tween::easeCubic(time, src.ystart, change.ystart, duration);
		yend = tween::easeCubic(time, src.yend, change.yend, duration);
		if (++time > duration) {
			return true;
		}
		return false;
	}
	
	void goToX(x_type start, x_type end) {
		if (start >= xmin) xstart = start;
		if (end <= xmax) xend = end;
		post();
	}
	
	void goToY(y_type start, y_type end) {
		if (start >= ymin) ystart = start;
		if (end <= ymax) yend = end;
		post();
	}
	
private:
	
	Window &parent;
	GLuint lists, nlists;
	x_type xstart, xend, xstep, xtickmajor, xease;
	
	const x_type xmin, xmax, xintervals;
	y_type ystart, yend, ystep, ytickmajor, yease;
	
	const y_type ymin, ymax, yintervals;
	const float tickSize, labelOffset, titleOffset, axisOffset, bgOverflow;
	const unsigned delay;
	float time;
	const float duration;
	
	GraphArea src, change;
	
	typedef std::list<plotter> plotlist;
	plotlist plots;
	
	std::queue<task> tasks;
	std::queue<GraphArea> dests;
	
	bool complete;
	
	
};

#endif
