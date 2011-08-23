#include "Window.hpp"

std::vector<Graph::xtype> g_x;
std::vector<Graph::ytype> g_y;

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
	ftlayout.SetFont(&font);
	
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
	
	graph = new Graph(*this);
	
	Graph::ytype y = 0;
	g_x.reserve(1000);
	g_y.reserve(1000);
	for (std::size_t i = 0; i < 1000; ++i) {
		g_x.push_back(std::rand() % 10000);
		y += double(std::rand()) / RAND_MAX * 0.1  - 0.05;
		g_y.push_back(y);
		//std::cout << g_x[g_x.size()-1] << ", " << g_y[g_y.size()-1] << std::endl;
	}
	std::sort(g_x.begin(), g_x.end());
	
	return true;
}

void Window::render() {
	// clear the screena nd depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// reset the view
	glLoadIdentity();
	
	glBegin(GL_QUADS);
	glColor3f(1, 0.9, 0.9); glVertex3f(0, 0, 0);
	glColor3f(0.9, 1, 0.9); glVertex3f(1, 0, 0);
	glColor3f(0.9, 0.9, 1); glVertex3f(1, 1, 0);
	glColor3f(1, 1, 1); glVertex3f(0, 1, 0);
	glEnd();
	
	/*
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
	*/
	
	// must call glColor before glRasterPos
	glColor3f(0.2, 0, 0);
	renderText("A journey of a thousand miles begins with a single step.", 0.5, 0.95);
	
	if (graph != NULL) {
		glTranslatef(0.15, 0.15, 0);
		glScalef(0.7, 0.7, 1);
		glColor3f(0.2, 0.2, 0.2);
		
		glCallList(graph->begin());
		graph->plot(g_x, g_y);
		
		glLoadIdentity();
		
	}
	
	
	SDL_GL_SwapBuffers();
}

void Window::onEvent(SDL_Event* event) {
	if (event->type == SDL_QUIT) {
		active = false;
	}
}

void Window::onLoop() {
	
}
