#include "Graph.hpp"

void Graph::init() {
	lists = glGenLists(nlists);
	
	char buf[20];
	
	glNewList(lists, GL_COMPILE);
	
	/// background
	{
		glColor3f(0.95, 0.95, 0.95);
		glBegin(GL_QUADS);
		glVertex2f(-bgOverflow, -bgOverflow);
		glVertex2f(-bgOverflow, height+bgOverflow);
		glVertex2f(width+bgOverflow, height+bgOverflow);
		glVertex2f(width+bgOverflow, -bgOverflow);
		glEnd();
	}
	
	/// grids
	{
		glColor3f(1, 1, 1);
		glLineWidth(2);
		
		for (x_type x = xbegin; x <= xend; x+=xtickmajor) {
			float xf = double(x-xbegin)/(xend-xbegin) * width;
			glBegin(GL_LINES);
			glVertex2f(xf, -bgOverflow);
			glVertex2f(xf, height+bgOverflow);
			glEnd();
		}
		
		for (y_type y = ybegin; y <= yend; y+=ytickmajor) {
			float yf = float(y-ybegin)/(yend-ybegin) * height;
			glBegin(GL_LINES);
			glVertex2f(-bgOverflow, yf);
			glVertex2f(width+bgOverflow, yf);
			glEnd();
		}
	}
	
	
	glColor3f(.3, .3, .3);
	
	/// x-axis
	{
		glTranslatef(0, -axisOffset, 0);
		
		// main line
		glLineWidth(2);
		glBegin(GL_LINES);
		glVertex2f(0, 0);
		glVertex2f(width, 0);
		glEnd();
		
		// title
		parent.renderText("Genomic coodinate", width/2, -titleOffset, 0, 18);

		glTranslatef(0, axisOffset, 0);
	}
	
	/// y-axis
	{
		glTranslatef(-axisOffset, 0, 0);
		
		// main line
		glLineWidth(2);
		glBegin(GL_LINES);
		glVertex2f(0, 0);
		glVertex2f(0, height);
		glEnd();
		
		// title
		parent.renderText("Relative copy number", -titleOffset, height/2, 90, 18);
		
		// ticks and labels
		for (y_type y = ybegin; y <= yend; y+=ytickmajor) {
			float yf = float(y-ybegin)/(yend-ybegin) * height;
			glBegin(GL_LINES);
			glVertex2f(tickSize, yf);
			glVertex2f(-tickSize, yf);
			glEnd();
			std::sprintf(buf, "%.1f", y);
			parent.renderText(buf, -(tickSize+labelOffset*0.5), yf, 0, 18, text::right);
		}
		
		glTranslatef(axisOffset, 0, 0);
	}
	
	/// reference line
	{
		glColor3f(0.8, 0.8, 0.8);
		glLineWidth(3);
		
		glBegin(GL_LINES);
		y_type zero_yf = -ybegin/float(yend-ybegin) * height;
		glVertex2f(-bgOverflow, zero_yf);
		glVertex2f(width+bgOverflow, zero_yf);
		glEnd();
	}
	
	glEndList();
}

void Graph::plot(const std::vector<x_type>& x, const std::vector<y_type>& y) {
	
	//scroll_x();
	
	easeto_x(3000, 5000);
	
	/// x axis
	if (skipflags ^ draw::xaxis) {
		char buf[20];
		
		glColor3f(.3, .3, .3);
		glLineWidth(2);
			
		glTranslatef(0, -axisOffset, 0);
			
		// ticks and labels
		for (x_type x = xbegin; x <= xend; x+=xtickmajor) {
			float xf = double(x-xbegin)/(xend-xbegin) * width;
			glBegin(GL_LINES);
			glVertex2f(xf, tickSize);
			glVertex2f(xf, -tickSize);
			glEnd();
			std::sprintf(buf, "%d", x);
			parent.renderText(buf, xf, -(tickSize+labelOffset), 0, 18);
		}
		
		glTranslatef(0, axisOffset, 0);
	}
	
	if (x.size() == 0) return;
	if (x.size() != y.size()) {
		throw std::runtime_error("x and y cannot differ in size");
	}
	
	float rangexf = float(xend-xbegin);
	float rangeyf = float(yend-ybegin);
	
	/// data lines
	{
		glLineWidth(2);
		glColor3f(.6, .6, .6);
		glPointSize(1);
		glBegin(GL_LINES);
		for (std::size_t i = 0; i < x.size()-1; ++i) {
			if (x[i] >= xbegin && x[i+1] <= xend) {
				glVertex2f(
					(x[i]-xbegin)/rangexf * width,
					(y[i]-ybegin)/rangeyf * height
				);
				glVertex2f(
					(x[i+1]-xbegin)/rangexf * width,
					(y[i+1]-ybegin)/rangeyf * height
				);
			}
		}
		glEnd();
	}
	
	/// data points
	{
		glPointSize(2);
		glBegin(GL_POINTS);
		for (std::size_t i = 0; i < x.size(); ++i) {
			if (x[i] >= xbegin && x[i] <= xend) {
				if (y[i] < -0.1) {
					glColor3f(0, 0, 1);
				} else if (y[i] > 0.1) {
					glColor3f(1, 0, 0);
				} else {
					glColor3f(0, 0, 0);
				}
				glVertex2f(
					(x[i]-xbegin)/rangexf * width,
					(y[i]-ybegin)/rangeyf * height
				);
			}
		}
		glEnd();
	}
	
}
	