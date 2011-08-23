#include "Graph.hpp"

void Graph::init() {
	lists = glGenLists(nlists);
	
	const float tickSize = 0.02;
	const float labelOffset = 0.03;
	const float titleOffset = 0.12;
	
	char buf[20];
	
	glNewList(lists, GL_COMPILE);
	
	/// x-axis
	
	// main line
	glBegin(GL_LINES);
	glVertex2f(0, 0);
	glVertex2f(1, 0);
	glEnd();
	
	// title
	parent.renderText("Genomic coodinate", 0.5, -titleOffset, 0, 18);

	// ticks and labels
	for (xtype x = xmin; x <= xmax; x+=xstep) {
		float xf = double(x)/xmax;
		glBegin(GL_LINES);
		glVertex2f(xf, 0);
		glVertex2f(xf, -tickSize);
		glEnd();
		std::sprintf(buf, "%d", x);
		parent.renderText(buf, xf, -(tickSize+labelOffset), 0, 18);
	}
	
	/// y-axis
	
	// main line
	glBegin(GL_LINES);
	glVertex2f(0, 0);
	glVertex2f(0, 1);
	glEnd();
	
	// title
	parent.renderText("Relative copy number", -titleOffset, 0.5, 90, 18);
	
	// ticks and labels
	for (ytype y = ymin; y <= ymax; y+=ystep) {
		float yf = float(y-ymin)/(ymax-ymin);
		glBegin(GL_LINES);
		glVertex2f(0, yf);
		glVertex2f(-tickSize, yf);
		glEnd();
		std::sprintf(buf, "%.1f", y);
		parent.renderText(buf, -(tickSize+labelOffset*0.5), yf, 0, 18, text::right);
	}
	
	/// zero reference line
	
	glBegin(GL_LINES);
	ytype zero_yf = -ymin/float(ymax-ymin);
	glVertex2f(0, zero_yf);
	glVertex2f(1, zero_yf);
	glEnd();
	
	
	glEndList();
}
