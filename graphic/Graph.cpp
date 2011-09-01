#include "Graph.hpp"

void Graph::init() {
	post();
	
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
		
		for (x_type x = xstart; x <= xend; x+=xtickmajor) {
			float xf = double(x-xstart)/(xend-xstart) * width;
			glBegin(GL_LINES);
			glVertex2f(xf, -bgOverflow);
			glVertex2f(xf, height+bgOverflow);
			glEnd();
		}
		
		for (y_type y = ystart; y <= yend; y+=ytickmajor) {
			float yf = float(y-ystart)/(yend-ystart) * height;
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
		
		glTranslatef(axisOffset, 0, 0);
	}
	
	
	glEndList();
	
	
	/*
	jumpTo(3000, 5000, ystart, yend);
	moveTo(1000, 2000, ystart, yend);
	easeTo(7000, 7500, ystart, yend);
	easeTo(0, 10000, ystart, yend);
	*/
	
}

void Graph::drawTicksX() {
	glTranslatef(0, -axisOffset, 0);
	
	glColor3f(.3, .3, .3);
	glLineWidth(2);
	
	// ticks and labels
	char buf[20];
	for (x_type x = xstart; x <= xend; x+=xtickmajor) {
		float xf = double(x-xstart)/(xend-xstart) * width;
		glBegin(GL_LINES);
		glVertex2f(xf, tickSize);
		glVertex2f(xf, -tickSize);
		glEnd();
		std::sprintf(buf, "%d", x);
		parent.renderText(buf, xf, -(tickSize+labelOffset), 0, 18);
	}
	glTranslatef(0, axisOffset, 0);
}

void Graph::drawTicksY() {
	glTranslatef(-axisOffset, 0, 0);
	
	glColor3f(.3, .3, .3);
	glLineWidth(2);
	
	// ticks and labels
	char buf[20];
	for (y_type y = ystart; y <= yend; y+=ytickmajor) {
		float yf = float(y-ystart)/(yend-ystart) * height;
		glBegin(GL_LINES);
		glVertex2f(tickSize, yf);
		glVertex2f(-tickSize, yf);
		glEnd();
		std::sprintf(buf, "%.1f", y);
		parent.renderText(buf, -(tickSize+labelOffset*0.5), yf, 0, 18, text::right);
	}
	
	glTranslatef(axisOffset, 0, 0);
}

void Graph::plot() {
	update();
	
	/// x axis
	if (draws.xaxis) {
		drawTicksX();
	}
	
	/// y axis
	if (draws.yaxis) {
		drawTicksY();
	}
	
	/// reference line
	if (draws.reference_line) {
		glColor3f(0.8, 0.8, 0.8);
		glLineWidth(3);
		
		glBegin(GL_LINES);
		y_type zero_yf = -ystart/float(yend-ystart) * height;
		glVertex2f(-bgOverflow, zero_yf);
		glVertex2f(width+bgOverflow, zero_yf);
		glEnd();
	}
	
	/// data
	plotlist::const_iterator it, end = plots.end();
	for (it = plots.begin(); it != end; ++it) {
		// call plotting function
		(*it)();
	}
	
}

void Graph::add(const std::vector<x_type>& x, const std::vector<y_type>& y) {
	if (x.size() == 0 || y.size() == 0) {
		throw std::runtime_error("x and y cannot be empty");
	}
	if (x.size() != y.size()) {
		throw std::runtime_error("x and y cannot differ in size");
	}
	
	// push pointer to member function onto plots queue, with parameters bound
	plots.push_back( boost::bind(&Graph::plot, this, x, y) );
}

void Graph::plot(const std::vector<x_type>& x, const std::vector<y_type>& y) {
	
	float rangexf = float(xend-xstart);
	float rangeyf = float(yend-ystart);
	
	/// data lines
	if (draws.data_lines) {
		glLineWidth(2);
		glColor3f(.6, .6, .6);
		glPointSize(1);
		glBegin(GL_LINES);
		for (std::size_t i = 0; i < x.size()-1; ++i) {
			if (x[i] >= xstart && x[i+1] <= xend) {
				glVertex2f(
					(x[i]-xstart)/rangexf * width,
					(y[i]-ystart)/rangeyf * height
				);
				glVertex2f(
					(x[i+1]-xstart)/rangexf * width,
					(y[i+1]-ystart)/rangeyf * height
				);
			}
		}
		glEnd();
	}
	
	/// data points
	if (draws.data_points) {
		glPointSize(2);
		glBegin(GL_POINTS);
		for (std::size_t i = 0; i < x.size(); ++i) {
			if (x[i] >= xstart && x[i] <= xend) {
				if (y[i] < -0.1) {
					glColor3f(0, 0, 1);
				} else if (y[i] > 0.1) {
					glColor3f(1, 0, 0);
				} else {
					glColor3f(0, 0, 0);
				}
				glVertex2f(
					(x[i]-xstart)/rangexf * width,
					(y[i]-ystart)/rangeyf * height
				);
			}
		}
		glEnd();
	}
	
}
	