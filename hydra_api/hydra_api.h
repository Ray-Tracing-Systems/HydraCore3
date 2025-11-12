// main header for HydraAPI 2.0
#pragma once
#include <cstdint>

#define HAPI_SEE_PUGIXML 
#ifdef HAPI_SEE_PUGIXML
  #include "pugixml.hpp"
#endif

#if defined(WIN32)
  #define HAPI
  //#ifdef HAPI_DLL
  //#define HAPI __declspec(dllexport) ///< mark all functions as 'extern "C"'; This is needed if you want to load DLL in dynamic;
  //#else
  //#define HAPI __declspec(dllimport) ///< mark all functions as 'extern "C"'; This is needed if you want to load DLL in dynamic;
  //#endif
#else
  //#if defined(__GNUC__)
  #define HAPI  ///< mark all functions as 'extern "C"'; This is needed if you want to load DLL dynamically;
#endif

struct HR2_SceneLibraryRef ///< main data storage, scene library, resource manager, navel of the earth. Usually this object is created in a single copy.
{              
  int32_t id = -1;
};

enum HR2_RES_STORAGE_TYPE { HR2_STORAGE_CPU  = 0, ///< force CPU implementation
                            HR2_STORAGE_GPU  = 1, ///< force GPU implementation, try not to store data on CPU when possible
                            HR2_STORAGE_DUAL = 2, ///< create both, duplicate data
};  

struct HR2_ReserveOpions
{
  int32_t maxGeoms     = 16384;
  int32_t maxImages    = 256;
  int32_t maxmaterials = 1024;
  int32_t maxLights    = 131126;
};

HAPI HR2_SceneLibraryRef hr2CreateLibrary(HR2_RES_STORAGE_TYPE a_type, HR2_ReserveOpions a_reserveOptions);
HAPI HR2_SceneLibraryRef hr2CreateLibraryFromFile(HR2_RES_STORAGE_TYPE a_type, HR2_ReserveOpions a_reserveOptions, const char* a_filename, bool a_async = false);

HAPI void hr2SaveSceneLibrary  (HR2_SceneLibraryRef, const char* a_filename, bool a_async = false);
HAPI void hr2DeleteSceneLibrary(HR2_SceneLibraryRef); ///< detele all
HAPI bool hr2SceneLibraryIsFinished(HR2_SceneLibraryRef a_cmbBuff); ///< check whether async scene load/save is completed; use this function within a wait-sleep loop when large scene is loaded/saved
                                                                                ///< in the first version async load/save is not planned for implementation, always return true

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct HR2_GeomRef     { int32_t id = -1; };
struct HR2_MaterialRef { int32_t id = -1; };
struct HR2_LightRef    { int32_t id = -1; };
struct HR2_TextureRef  { int32_t id = -1; };
struct HR2_SpectrumRef { int32_t id = -1; };
struct HR2_CameraRef   { int32_t id = -1; };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum HR2_CMD_TYPE { HR2_CMDBUF_APPEND = 0, ///<! create new bjects
                    HR2_CMDBUF_UPDATE = 1, ///<! update existing objects
                    HR2_CMDBUF_DUAL   = 2, ///<! simultaniously create and update
                    HR2_CMDBUF_UNDEFINED = 0xFFFFFFFF,
}; 

struct HR2_CommandBuffer ///<! use this object to add new data to scene library
{
  int32_t id = -1;
  HR2_CMD_TYPE type = HR2_CMDBUF_UNDEFINED;
};

HAPI HR2_CommandBuffer hr2CreateCommandBuffer(HR2_SceneLibraryRef a_scnLib, HR2_CMD_TYPE a_type);
HAPI void              hr2CommitCommandBuffer(HR2_CommandBuffer a_cmbBuff, bool a_async = false); ///<! Commit and immediately delete it

HAPI bool              hr2CommandBufferIsFinished(HR2_CommandBuffer a_cmbBuff); ///<! check wherther async commit is completed; use this function within a wait-sleep loop when large scene is loaded/added
                                                                                   ///<! in the first version async commit is not planned for implementation, always return true

//// Create new objects 

HAPI HR2_MaterialRef hr2CreateMaterial(HR2_CommandBuffer a_cmdBuff, const char* a_matName = "");
HAPI HR2_LightRef    hr2CreateLight   (HR2_CommandBuffer a_cmdBuff, const char* a_lgtName = "");
HAPI HR2_CameraRef   hr2CreateCamera  (HR2_CommandBuffer a_cmdBuff, const char* a_camName = "");

struct HR2_MeshInput
{
  float*   vPosPtr     = nullptr; 
  uint32_t vPosStride  = 4;
  float*   vNormPtr    = nullptr; 
  uint32_t vNormStride = 4;
  float*   vTangPtr    = nullptr; 
  uint32_t vTangStride = 4;
  float*   vTexCoordPtr    = nullptr; 
  uint32_t vTexCoordStride = 2;

  uint32_t* indicesPtr = nullptr;
  uint32_t  indicesNum = 0;
  //uint32_t  indicesStride = 3; ///<! 3 for triangle mesh, 4 for quad mesh, quads are not supported yet

  uint32_t* matIdPtr = nullptr;
  uint32_t  matIdAll = 0;
  uint32_t  matIdNum = 1; ///<! if 1, set whole mesh with single material, read matIdAll; else read material indices from matIdPtr
};

HAPI HR2_GeomRef    hr2CreateMeshFromData(HR2_CommandBuffer a_cmdBuff, const char* a_meshName, HR2_MeshInput a_input);
HAPI HR2_GeomRef    hr2CreateGeomFromFile(HR2_CommandBuffer a_cmdBuff, const char* a_filename);
HAPI HR2_TextureRef hr2CreateTextureFromFile(HR2_CommandBuffer a_cmdBuff, const char* a_filename);

//// Update existing objects

HAPI HR2_GeomRef    hr2GetGeom   (HR2_CommandBuffer a_cmdBuff, const int32_t a_id); ///<! Geometry and texture update will not be implemented in first version
HAPI HR2_TextureRef hr2GetTexture(HR2_CommandBuffer a_cmdBuff, const int32_t a_id); ///<! Geometry and texture update will not be implemented in first version
 
#ifdef HAPI_SEE_PUGIXML
HAPI pugi::xml_node hr2MaterialParamNode(HR2_MaterialRef a_mat);
HAPI pugi::xml_node hr2LightParamNode(HR2_LightRef a_lgt);
HAPI pugi::xml_node hr2CameraParamNode(HR2_CameraRef a_cam);
#endif