#include "mesh_load_obj.h"

#include <cstdio>
#include <string>
#include <unordered_map>

#ifdef DISABLE_OBJ_LOADER
SimpleMesh LoadMeshFromObj(const char* a_fileName, bool aVerbose, const char* a_mtlBaseDir)
{
  printf("[LoadMeshFromObj::ERROR] OBJ loader is disabled\n");
  return SimpleMesh();
}
#else

#define TINYOBJLOADER_IMPLEMENTATION
#include "3rd_party/tiny_obj_loader.h"

struct TinyObjIndexEqual
{
  bool operator()(const tinyobj::index_t &lhs, const tinyobj::index_t &rhs) const
  {
    return lhs.vertex_index == rhs.vertex_index &&
           lhs.normal_index == rhs.normal_index &&
           lhs.texcoord_index == rhs.texcoord_index;
  }
};

struct TinyObjIndexHasher
{
  size_t operator()(const tinyobj::index_t &index) const
  {
    return ((std::hash<int>()(index.vertex_index) ^
            (std::hash<int>()(index.normal_index) << 1)) >> 1) ^
            (std::hash<int>()(index.texcoord_index) << 1);
  }
};

namespace cmesh4
{
  SimpleMesh LoadMeshFromObj(const char* a_fileName, bool aVerbose, const char* a_mtlBaseDir)
  {
    SimpleMesh mesh;

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    
    std::string warn;
    std::string err;

    bool loading_result = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, a_fileName, a_mtlBaseDir);

    if (!loading_result)
    {
      printf("[LoadMeshFromObj::ERROR] Failed to load obj file: %s\n", err.c_str());
      return mesh;
    }

    if (aVerbose)
    {
      if (warn.empty())
        printf("[LoadMeshFromObj::INFO] Loaded obj file: %s\n", a_fileName);
      else
        printf("[LoadMeshFromObj::WARNING] Loaded obj file %s with warnings: %s\n", a_fileName, warn.c_str());
    }

    const LiteMath::float4 default_norm = float4(0, 0, 1, 0);
    const LiteMath::float4 default_tangent = float4(1, 0, 0, 0);
    const LiteMath::float2 default_texcoord = float2(0, 0);

    std::unordered_map<tinyobj::index_t, uint32_t, TinyObjIndexHasher, TinyObjIndexEqual> uniqueVertIndices = {};

    uint32_t numIndices = 0;
    for (const auto& shape : shapes)
      numIndices += shape.mesh.indices.size();

    mesh.vPos4f.reserve(attrib.vertices.size() / 3);
    mesh.vNorm4f.reserve(attrib.vertices.size() / 3);
    mesh.vTang4f.reserve(attrib.vertices.size() / 3);
    mesh.vTexCoord2f.reserve(attrib.vertices.size() / 3);
    mesh.indices.reserve(numIndices);

    for (const auto& shape : shapes)
    {
      mesh.matIndices.insert(std::end(mesh.matIndices), std::begin(shape.mesh.material_ids), std::end(shape.mesh.material_ids));

      for (const auto& index : shape.mesh.indices)
      {
        uint32_t my_index = 0;
        auto it = uniqueVertIndices.find(index);
        if (it != uniqueVertIndices.end())
        {
          my_index = it->second;
        }
        else
        {
          uniqueVertIndices.insert({index, static_cast<uint32_t>(mesh.vPos4f.size())});
          my_index = static_cast<uint32_t>(mesh.vPos4f.size());

          assert(index.vertex_index >= 0 && index.vertex_index < attrib.vertices.size() / 3);
          mesh.vPos4f.push_back({attrib.vertices[3 * index.vertex_index + 0],
                                 attrib.vertices[3 * index.vertex_index + 1],
                                 attrib.vertices[3 * index.vertex_index + 2],
                                 1.0f});
          if(index.normal_index >= 0)
          {
            mesh.vNorm4f.push_back({attrib.normals[3 * index.normal_index + 0],
                                    attrib.normals[3 * index.normal_index + 1],
                                    attrib.normals[3 * index.normal_index + 2],
                                    0.0f});
          }
          else
          {
            mesh.vNorm4f.push_back(default_norm);
          }
          if(index.texcoord_index >= 0)
          {
            mesh.vTexCoord2f.push_back({attrib.texcoords[2 * index.texcoord_index + 0],
                                        attrib.texcoords[2 * index.texcoord_index + 1]});
          }
          else
          {
            mesh.vTexCoord2f.push_back(default_texcoord);
          }
          mesh.vTang4f.push_back(default_tangent);
        }

        mesh.indices.push_back(my_index);
      }
    }

    // fix material id
    for (unsigned &mid : mesh.matIndices) 
    { 
      if(mid == uint32_t(-1)) 
        mid = 0;
    }

    if (aVerbose)
    {
      printf("[LoadMeshFromObj::INFO] Loaded obj file %s with %d vertices and %d indices\n", 
              a_fileName, (unsigned)mesh.vPos4f.size(), (unsigned)mesh.indices.size());
    }

    return mesh;
  }
}
#endif