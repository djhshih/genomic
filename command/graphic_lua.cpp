/*
** Lua binding: graphic
** Generated automatically by tolua++-1.0.92 on Wed Aug 31 21:24:10 2011.
*/

#ifndef __cplusplus
#include "stdlib.h"
#endif
#include "string.h"

#include "tolua++.h"

/* Exported function */
TOLUA_API int  tolua_graphic_open (lua_State* tolua_S);

#include "graphic.hpp"

/* function to register type */
static void tolua_reg_types (lua_State* tolua_S)
{
}

/* function: jumpTo */
#ifndef TOLUA_DISABLE_tolua_graphic_jumpTo00
static int tolua_graphic_jumpTo00(lua_State* tolua_S)
{
#ifndef TOLUA_RELEASE
 tolua_Error tolua_err;
 if (
     !tolua_isnumber(tolua_S,1,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,2,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,3,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,4,0,&tolua_err) ||
     !tolua_isnoobj(tolua_S,5,&tolua_err)
 )
  goto tolua_lerror;
 else
#endif
 {
   long xstartDest = ((  long)  tolua_tonumber(tolua_S,1,0));
   long xendDest = ((  long)  tolua_tonumber(tolua_S,2,0));
   float ystartDest = ((  float)  tolua_tonumber(tolua_S,3,0));
   float yendDest = ((  float)  tolua_tonumber(tolua_S,4,0));
  {
   jumpTo(xstartDest,xendDest,ystartDest,yendDest);
  }
 }
 return 0;
#ifndef TOLUA_RELEASE
 tolua_lerror:
 tolua_error(tolua_S,"#ferror in function 'jumpTo'.",&tolua_err);
 return 0;
#endif
}
#endif //#ifndef TOLUA_DISABLE

/* function: moveTo */
#ifndef TOLUA_DISABLE_tolua_graphic_moveTo00
static int tolua_graphic_moveTo00(lua_State* tolua_S)
{
#ifndef TOLUA_RELEASE
 tolua_Error tolua_err;
 if (
     !tolua_isnumber(tolua_S,1,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,2,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,3,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,4,0,&tolua_err) ||
     !tolua_isnoobj(tolua_S,5,&tolua_err)
 )
  goto tolua_lerror;
 else
#endif
 {
   long xstartDest = ((  long)  tolua_tonumber(tolua_S,1,0));
   long xendDest = ((  long)  tolua_tonumber(tolua_S,2,0));
   float ystartDest = ((  float)  tolua_tonumber(tolua_S,3,0));
   float yendDest = ((  float)  tolua_tonumber(tolua_S,4,0));
  {
   moveTo(xstartDest,xendDest,ystartDest,yendDest);
  }
 }
 return 0;
#ifndef TOLUA_RELEASE
 tolua_lerror:
 tolua_error(tolua_S,"#ferror in function 'moveTo'.",&tolua_err);
 return 0;
#endif
}
#endif //#ifndef TOLUA_DISABLE

/* function: easeTo */
#ifndef TOLUA_DISABLE_tolua_graphic_easeTo00
static int tolua_graphic_easeTo00(lua_State* tolua_S)
{
#ifndef TOLUA_RELEASE
 tolua_Error tolua_err;
 if (
     !tolua_isnumber(tolua_S,1,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,2,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,3,0,&tolua_err) ||
     !tolua_isnumber(tolua_S,4,0,&tolua_err) ||
     !tolua_isnoobj(tolua_S,5,&tolua_err)
 )
  goto tolua_lerror;
 else
#endif
 {
   long xstartDest = ((  long)  tolua_tonumber(tolua_S,1,0));
   long xendDest = ((  long)  tolua_tonumber(tolua_S,2,0));
   float ystartDest = ((  float)  tolua_tonumber(tolua_S,3,0));
   float yendDest = ((  float)  tolua_tonumber(tolua_S,4,0));
  {
   easeTo(xstartDest,xendDest,ystartDest,yendDest);
  }
 }
 return 0;
#ifndef TOLUA_RELEASE
 tolua_lerror:
 tolua_error(tolua_S,"#ferror in function 'easeTo'.",&tolua_err);
 return 0;
#endif
}
#endif //#ifndef TOLUA_DISABLE

/* Open function */
TOLUA_API int tolua_graphic_open (lua_State* tolua_S)
{
 tolua_open(tolua_S);
 tolua_reg_types(tolua_S);
 tolua_module(tolua_S,NULL,0);
 tolua_beginmodule(tolua_S,NULL);
  tolua_function(tolua_S,"jumpTo",tolua_graphic_jumpTo00);
  tolua_function(tolua_S,"moveTo",tolua_graphic_moveTo00);
  tolua_function(tolua_S,"easeTo",tolua_graphic_easeTo00);
 tolua_endmodule(tolua_S);
 return 1;
}


#if defined(LUA_VERSION_NUM) && LUA_VERSION_NUM >= 501
 TOLUA_API int luaopen_graphic (lua_State* tolua_S) {
 return tolua_graphic_open(tolua_S);
};
#endif

