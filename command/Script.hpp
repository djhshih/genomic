#ifndef command_SCRIPT_h
#define command_SCRIPT_h

#include <cstdlib>
#include <cstring>
#include <stdexcept>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

// Scripting engine


using namespace std;

class Script
{
public:
	
	int exec() {
		// create state
		lua_State *L = lua_open();
		if (L == NULL) {
			throw std::runtime_error("Failed to initialize Lua interpreter.");
		}

		// stop collector during initialization
		lua_gc(L, LUA_GCSTOP, 0);
		// open libraries
		luaL_openlibs(L);
		// restart collector
		lua_gc(L, LUA_GCRESTART, 0);

		char prompt[] = "> ";

		const unsigned maxlinelen = 256;
		char line[maxlinelen];
		int error;
		
		printf(prompt);
		while (fgets(line, maxlinelen-1, stdin) != NULL) {
			error = luaL_loadbuffer(L, line, strlen(line), "") ||
				lua_pcall(L, 0, 0, 0);
			if (error) {
				// print error message from stack
				fprintf(stderr, "%s\n", lua_tostring(L, -1));
				// pop message from stack
				lua_pop(L, 1);
			}
			printf(prompt);
		}

		// clean up
		lua_close(L);
	}
	
private:
	
	
};

#endif
