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
		bool background;
		PropertiesDraw()
		: xaxis(true), yaxis(true), background(true) {}
	};
	
	struct Area {
		x_type xbegin, xend;
		y_type ybegin, yend;
		Area(x_type xb, x_type xe, y_type yb, y_type ye)
		: xbegin(xb), xend(xe), ybegin(yb), yend(ye) {}
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
	  xbegin(0), xend(10000), xstep(1), xintervals(10),
	  ybegin(-2), yend(4), ystep(0.5), yintervals(6),
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
		xbegin += xstep;
		xend += xstep;
		if (xbegin < xmin) xbegin = xmin;
		if (xend > xmax) xend = xmax;
		SDL_Delay(delay);
	}
	
	void scrollY() {
		ybegin += ystep;
		yend += ystep;
		if (ybegin < ymin) ybegin = ymin;
		if (yend > ymax) yend = ymax;
		SDL_Delay(delay);
	}
	
	void addTask(task t, x_type xbeginDest, x_type xendDest, y_type ybeginDest, y_type yendDest) {
		dests.push(Area(xbeginDest, xendDest, ybeginDest, yendDest));
		tasks.push(t);
	}
	
	void jumpTo(x_type xbeginDest, x_type xendDest, y_type ybeginDest, y_type yendDest) {
		addTask(&Graph::jumpToTask, xbeginDest, xendDest, ybeginDest, yendDest);
	}
	
	bool jumpToTask() {
		xbegin += xbegin_change;
		xend += xend_change;
		ybegin += ybegin_change;
		yend += yend_change;
		if (xbegin < xmin) xbegin = xmin;
		if (ybegin < ymin) ybegin = ymin;
		if (xend > xmax) xend = xmax;
		if (yend > ymax) yend = ymax;
		return true;
	}
	
	void moveTo(x_type xbeginDest, x_type xendDest, y_type ybeginDest, y_type yendDest) {
		addTask(&Graph::moveToTask, xbeginDest, xendDest, ybeginDest, yendDest);
	}
	
	bool moveToTask() {
		xbegin = tween::linear(time, xbegin_src, xbegin_change, duration);
		xend = tween::linear(time, xend_src, xend_change, duration);
		ybegin = tween::linear(time, ybegin_src, ybegin_change, duration);
		yend = tween::linear(time, yend_src, yend_change, duration);
		if (++time > duration) {
			return true;
		}
		return false;
	}
	
	void easeTo(x_type xbeginDest, x_type xendDest, y_type ybeginDest, y_type yendDest) {
		addTask(&Graph::easeToTask, xbeginDest, xendDest, ybeginDest, yendDest);
	}
	
	bool easeToTask() {
		xbegin = tween::easeCubic(time, xbegin_src, xbegin_change, duration);
		xend = tween::easeCubic(time, xend_src, xend_change, duration);
		ybegin = tween::easeCubic(time, ybegin_src, ybegin_change, duration);
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
		taskqueue::iterator it = tasks.begin();
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
			xbegin_change = dest.xbegin - xbegin;
			xend_change = dest.xend - xend;
			ybegin_change = dest.ybegin - ybegin;
			yend_change = dest.yend - yend;
			
			xbegin_src = xbegin;
			xend_src = xend;
			ybegin_src = ybegin;
			yend_src = yend;
			
			time = 0;
			complete = false;
		}
	}
	
	void post() {
		xtickmajor = (xend - xbegin) / xintervals;
		ytickmajor = (yend - ybegin) / yintervals;
		complete = true;
	}
	
	/*
	void easeto_x(x_type begin, x_type end) {
		if (begin < xmin) begin = xmin;
		if (end > xmax) end = xmax;
		
		xbegin += (begin - xbegin) / xease;
		xend += (end - xend) / xease;
		xtickmajor = (xend - xbegin) / xintervals;
		
		if (abs(xbegin - begin) < xease) {
			xbegin = begin;
			xend = end;
			draws.xaxis = true;
		} else {
			draws.xaxis = false;
		}
		
		SDL_Delay(delay);
	}
	*/
	
	void plot(const std::vector<x_type>& x, const std::vector<y_type>& y);
	
	GLuint begin() {
		return lists;
	}
	
	GLuint end() {
		return lists + nlists;
	}
	
private:
	
	void init();
	
	void jumpToX(x_type begin, x_type end) {
		if (begin >= xmin) xbegin = begin;
		if (end <= xmax) xend = end;
		post();
	}
	
	void jumpToY(y_type begin, y_type end) {
		if (begin >= ymin) ybegin = begin;
		if (end <= ymax) yend = end;
		post();
	}
	
private:
	
	Window &parent;
	GLuint lists, nlists;
	x_type xbegin, xend, xstep, xtickmajor, xease;
	x_type xbegin_src, xbegin_change, xend_src, xend_change;
	const x_type xmin, xmax, xintervals;
	y_type ybegin, yend, ystep, ytickmajor, yease;
	y_type ybegin_src, ybegin_change, yend_src, yend_change;
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
