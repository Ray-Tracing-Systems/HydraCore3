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

struct HAPI_SceneLibrary ///< main data storage, scene library, resource manager, navel of the earth. Usually this object is created in a single copy.
{              
  int32_t id = -1;
};

enum HAPI_RES_STORAGE_TYPE { HAPI_STORAGE_CPU  = 0, ///< force CPU implementation
                             HAPI_STORAGE_GPU  = 1, ///< force GPU implementation, try not to store data on CPU when possible
                             HAPI_STORAGE_DUAL = 2, ///< create both, duplicate data
};  

struct HAPI_ReserveOpions
{
  int32_t maxGeoms     = 16384;
  int32_t maxImages    = 256;
  int32_t maxmaterials = 1024;
  int32_t maxLights    = 131126;
};

HAPI HAPI_SceneLibrary hapiCreateLibraryEmpty   (HAPI_RES_STORAGE_TYPE a_type, HAPI_ReserveOpions a_reserveOptions);
HAPI HAPI_SceneLibrary hapiCreateLibraryFromFile(HAPI_RES_STORAGE_TYPE a_type, HAPI_ReserveOpions a_reserveOptions, const char* a_filename, bool a_async = false);

HAPI void              hapiSaveSceneLibrary  (HAPI_SceneLibrary, const char* a_filename, bool a_async = false);
HAPI void              hapiDeleteSceneLibrary(HAPI_SceneLibrary); ///< detele all

HAPI bool              hapiSceneLibraryIsFinished(HAPI_SceneLibrary a_cmbBuff); ///< check whether async scene load/save is completed; use this function within a wait-sleep loop when large scene is loaded/saved
                                                                                ///< in the first version async load/save is not planned for implementation, always return true

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct HAPI_Geom     { int32_t id = -1; };
struct HAPI_Material { int32_t id = -1; };
struct HAPI_Light    { int32_t id = -1; };
struct HAPI_Texture  { int32_t id = -1; };
struct HAPI_Spectrum { int32_t id = -1; };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum HAPI_CMD_TYPE { HAPI_CMD_APPEND = 0, ///<! create new bjects
                     HAPI_CMD_UPDATE = 1, ///<! update existing objects
                     HAPI_CMD_DUAL   = 2, ///<! simultaniously create and update
                     HAPI_CMD_UNDEFINED = 0xFFFFFFFF,
}; 

struct HAPI_CommandBuffer ///<! use this object to add new data to scene library
{
  int32_t id = -1;
  HAPI_CMD_TYPE type = HAPI_CMD_UNDEFINED;
};

HAPI HAPI_CommandBuffer hapiCreateCommandBuffer(HAPI_SceneLibrary a_scnLib, HAPI_CMD_TYPE a_type);
HAPI void               hapiCommitCommandBuffer(HAPI_CommandBuffer a_cmbBuff, bool a_async = false); ///<! Commit and immediately delete it

HAPI bool               hapiCommandBufferIsFinished(HAPI_CommandBuffer a_cmbBuff); ///<! check wherther async commit is completed; use this function within a wait-sleep loop when large scene is loaded/added
                                                                                   ///<! in the first version async commit is not planned for implementation, always return true

//// Create new objects 

HAPI HAPI_Material hapiCreateMaterialEmpty  (HAPI_CommandBuffer a_cmdBuff, const char* a_matName = "");


HAPI HAPI_Geom     hapiCreateMeshEmpty      (HAPI_CommandBuffer a_cmdBuff, const char* a_meshName);
HAPI HAPI_Geom     hapiCreateGeomFromFile   (HAPI_CommandBuffer a_cmdBuff, const char* a_filename);
HAPI HAPI_Texture  hapiCreateTextureFromFile(HAPI_CommandBuffer a_cmdBuff, const char* a_filename);

//// Update existing objects

HAPI HAPI_Geom    hapiGetGeom   (HAPI_CommandBuffer a_cmdBuff, const int32_t a_id); ///<! Geometry and texture update will not be implemented in first version
HAPI HAPI_Texture hapiGetTexture(HAPI_CommandBuffer a_cmdBuff, const int32_t a_id); ///<! Geometry and texture update will not be implemented in first version
 
#ifdef HAPI_SEE_PUGIXML
HAPI pugi::xml_node hapiMaterialParamNode(HAPI_Material a_mat);
#endif