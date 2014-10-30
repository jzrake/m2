#ifndef LUA_STRUCT_H
#define LUA_STRUCT_H
#include <lua.h>

enum {
  LSTRUCT_ORPHAN  = -99119911,
  LSTRUCT_CHAR    =  11000001,
  LSTRUCT_INT     =  12000001,
  LSTRUCT_ULONG   =  13000001,
  LSTRUCT_DOUBLE  =  14000001,
  LSTRUCT_STRING  =  15000001,
  LSTRUCT_OBJECT  =  16000001,
  LSTRUCT_STRUCT  =  17000001,
  LSTRUCT_POINTER =  18000001,
} ;


typedef struct {
  const char *method_name;
  lua_CFunction method;
} lua_struct_method_t;


typedef struct {
  const char *member_name;
  ptrdiff_t offset;
  int data_type;
  const char *type_name; /* if type is LSTRUCT_STRUCT */
  unsigned string_len; /* if type is LSTRUCT_STRING */
} lua_struct_member_t;


typedef struct {
  const char *type_name;
  size_t alloc_size;
  lua_struct_member_t *members;
  lua_struct_method_t *methods;
  lua_CFunction __tostring;
  lua_CFunction __len;
  lua_CFunction __call;
  lua_CFunction __init;
  lua_CFunction __gc;
} lua_struct_t;


int luaopen_struct(lua_State *L);
int lua_struct_register(lua_State *L, lua_struct_t type);
int lua_struct_pushparent(lua_State *L, int n);
int lua_struct_pushmember(lua_State *L, int n, const char *member_name);
void *lua_struct_new(lua_State *L, const char *tname);
void *lua_struct_newref(lua_State *L, const char *tname, void *val, int parent);
void *lua_struct_checkstruct(lua_State *L, int n, const char *tname);
lua_struct_t lua_struct_newtype(lua_State *L);

#endif /* LUA_STRUCT_H */
