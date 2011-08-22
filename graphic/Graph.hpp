#ifndef graphic_Graph_h
#define

#include <GL/gl.h>

class Graph
{
public:
	float width, height;
	
public:
	
	Graph() : nlists(1) {
		lists = glGenLists(nlists);
		
		glNewList(lists, GL_COMPILE);
		
		// x-axis
		glBegin(GL_LINES);
		glVertex2f(0, 0);
		glVertex2f(1, 0);
		glEnd();
		
		// y-axis
		glBegin(GL_LINES);
		glVertex2f(0, 0);
		glVertex2f(0, 1);
		glEnd();
		
		
		
		
		glEndList();
		
	}
	
	GLuint begin() {
		return lists;
	}
	GLuint end() {
		return lists + nlists;
	}
private:
	GLuint lists, nlists;
	
};

#endif