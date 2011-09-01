
#include "ScriptGraphic.hpp"
//#include "../graphic/WindowThreaded.hpp"
#include "../graphic/Window.hpp"
#include "graphic.hpp"


/*

/// Method 1

class ScriptableWindow : public Window {
	Script s;
	// note: Window::onLoop must be virtual
	void onLoop() {
		s.prompt();
	}
};

int main(int argc, char *argv[]) {
	ScriptableWindow sw;
	return sw.exec();
}
*/


/*

/// Method 2

struct ScriptCall : public Callable {
	Script s;
	void operator()() {
		s.prompt();
	}
};

int main(int argc, char *argv[]) {
	ScriptCall call;
	// pass functor to Window
	Window w(&call);
	return w.exec();
}
*/

/*
Script script;

int callback(void *) {
	script.exec();
}

int main(int argc, char *argv[]) {
	WindowThreaded w(&callback);
	return w.exec();
}
*/


ScriptGraphic script(window);

int callback(void *) {
	easeTo(2000, 6000, -2, 4);
	script.exec();
}

int main(int argc, char *argv[]) {
	window.addThread(&callback);
	return window.exec();
}
