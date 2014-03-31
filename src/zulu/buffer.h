#ifndef LUA_BUFFER_H
#define LUA_BUFFER_H

#include "lua.h"

enum { 
  BUFFER_TYPE_CHAR=1,
  BUFFER_TYPE_INT=2,
  BUFFER_TYPE_LONG=3,
  BUFFER_TYPE_FLOAT=4,
  BUFFER_TYPE_DOUBLE=5,
} ;

void *buf_new_buffer(lua_State *L, const void *p, size_t size);

#endif /* LUA_BUFFER_H */
