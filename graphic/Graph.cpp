#include "Graph.hpp"

void Graph::init() {
	lists = glGenLists(nlists);
	
	const float tickSize = 0.1;
	
	glNewList(lists, GL_COMPILE);
	
	// x-axis
	
	// main line
	glBegin(GL_LINES);
	glVertex2f(0, 0);
	glVertex2f(1, 0);
	glEnd();
	
	//glColor3f(1, 1, 1);
	glRasterPos2f(0.5, -0.1);
	parent.renderText("X-axis");

	// ticks
	for (long x = xmin; x <= xmax; x+=xstep) {
		glBegin(GL_LINES);
		float xf = double(x)/xmax;
		glVertex2f(xf, 0);
		glVertex2f(xf, -tickSize);
		glEnd();
	}
	
	// y-axis
	
	// main line
	glBegin(GL_LINES);
	glVertex2f(0, rescale(ymin));
	glVertex2f(0, rescale(ymax));
	glEnd();
	
	// ticks
	for (float y = ymin; y <= ymax; y+=ystep) {
		glBegin(GL_LINES);
		float yf = rescale(y);
		glVertex2f(0, yf);
		glVertex2f(-tickSize, yf);
		glEnd();
	}
	
	glEndList();
}
