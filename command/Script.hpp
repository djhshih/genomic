#ifndef command_SCRIPT_h
#define command_SCRIPT_h

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
using namespace std;

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

// Scripting engine


class Script
{
public:
	
	Script()
	: L(NULL), promptStr("> ") {
		init();
	}
	
	~Script() {
		// clean up
		if (L != NULL) lua_close(L);
	}
	
	void init() {
		// create state
		L = lua_open();
		if (L == NULL) {
			throw std::runtime_error("Failed to initialize Lua interpreter.");
		}

		// stop collector during initialization
		lua_gc(L, LUA_GCSTOP, 0);
		// open libraries
		luaL_openlibs(L);
		// restart collector
		lua_gc(L, LUA_GCRESTART, 0);
		
	}
	
	int exec() {
		while (prompt());
	}
	
	bool prompt() {
		int error;
		printf(promptStr.c_str());
		if (fgets(line, maxlinelen-1, stdin) != NULL) {
			error = luaL_loadbuffer(L, line, strlen(line), "") ||
				lua_pcall(L, 0, 0, 0);
			if (error) {
				// print error message from stack
				fprintf(stderr, "%s\n", lua_tostring(L, -1));
				// pop message from stack
				lua_pop(L, 1);
			}
			return true;
		} else {
			return false;
		}
	}
	
private:
	
	lua_State *L;
	
	string promptStr;
	static const unsigned maxlinelen = 256;
	char line[maxlinelen];
	
};

#endif
