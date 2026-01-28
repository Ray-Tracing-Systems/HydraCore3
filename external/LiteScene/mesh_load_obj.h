#ifndef LITESCENE_MESH_LOAD_OBJ_H_
#define LITESCENE_MESH_LOAD_OBJ_H_
#include "cmesh4.h"

namespace cmesh4
{
  SimpleMesh LoadMeshFromObj(const char* a_fileName, bool aVerbose = false, const char* a_mtlBaseDir = nullptr);
}

#endif