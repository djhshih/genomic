#include "Window.hpp"

int Window::exec() {
	if (!init()) return -1;
	
	SDL_Event event;
	while (active) {
		while (SDL_PollEvent(&event)) {
			onEvent(&event);
		}
		onLoop();
		render();
	}
	
	return 0;
}

bool Window::init() {
	// initialize SDL
	if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
			return false;
	}
	
	// set video mode
	surface = SDL_SetVideoMode(width, height, color, SDL_HWSURFACE | SDL_GL_DOUBLEBUFFER | SDL_OPENGL);
	if (surface == NULL) {
		return false;
	}
	
	// set title bar
	SDL_WM_SetCaption("Window", NULL);
	
	if (font.Error()) {
		std::cerr << "Failed to load font" << std::endl;
	}
	
	// set background colour
	glClearColor(0, 0, 0, 0);
	
	glViewport(0, 0, width, height);
	
	// configure the porject matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// left, right, bottom, top, near, far
	// set (0, 0) at bottom left
	glOrtho(0, 1, 0, 1, -1, 1);
	
	// configure the model view matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	glEnable(GL_TEXTURE_2D);
	
	return true;
}

void Window::render() {
	// clear the screena nd depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// reset the view
	glLoadIdentity();
	
	glBegin(GL_QUADS);
	glColor3f(1, 0, 0); glVertex3f(0, 0, 0);
	glColor3f(0, 1, 0); glVertex3f(1, 0, 0);
	glColor3f(0, 0, 1); glVertex3f(1, 1, 0);
	glColor3f(1, 1, 1); glVertex3f(0, 1, 0);
	glEnd();
	
	const int n = 10;
	
	GLuint list = glGenLists(2);
	
	glNewList(list, GL_COMPILE);
	glBegin(GL_LINES);
	glVertex2f(0, 0); glVertex2f(1, 0);
	glEnd();
	glEndList();
	
	glNewList(list+1, GL_COMPILE);
	glBegin(GL_LINES);
	glVertex2f(0, 0.01); glVertex2f(1, 0.01);
	glEnd();
	glEndList();
	
	for (int i = 0; i < n; ++i) {
		glTranslatef(0, i*0.05+0.1, 0);
		glColor3f(1, 1, 0);
		glCallList(list);
		glColor3f(0, 1, 1);
		glCallList(list+1);
		glLoadIdentity();
	}
	
	// must call glColor before glRasterPos
	glColor3f(0, 0, 0);
	glRasterPos2f(0.1, 0.9);
	font.FaceSize(12);
	font.Render("A journey of a thousand miles begins with a single step.");
		
	SDL_GL_SwapBuffers();
}

void Window::onEvent(SDL_Event* event) {
	if (event->type == SDL_QUIT) {
		active = false;
	}
}

void Window::onLoop() {
	
}