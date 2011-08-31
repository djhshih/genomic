
#include "Script.hpp"
#include "../graphic/WindowThreaded.hpp"


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

Script script;

int callback(void *) {
	script.exec();
}

int main(int argc, char *argv[]) {
	WindowThreaded w(&callback);
	return w.exec();
}
