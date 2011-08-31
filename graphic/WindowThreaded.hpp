#ifndef graphic_WindowThreaded_h
#define graphic_WindowThreaded_h

#include <stdexcept>

#include <SDL/SDL_thread.h>

#include "Window.hpp"

class WindowThreaded : public Window
{
	typedef int (*thread_callback)(void *);
	thread_callback callback;
	
public:
	
	WindowThreaded()
	: callback(NULL) {
		init();
	}
	
	WindowThreaded(thread_callback _callback)
	: callback(_callback) {
		init();
	}
	
private:
	
	void init() {
		/// Setup thread
		if (callback != NULL) {
			SDL_Thread *thread = SDL_CreateThread(callback, NULL);
			if (thread == NULL) {
				throw std::runtime_error("Failed to create thread.");
			}
		}
	}
};

#endif
