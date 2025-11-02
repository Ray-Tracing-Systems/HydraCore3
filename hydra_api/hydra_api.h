// main header for HydraAPI 2.0
#pragma once
#include <cstdint>

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

enum HAPI_RES_MANAGER_TYPE { HAPI_RES_MGR_CPU     = 0, ///< force CPU implementation
                             HAPI_RES_MGR_GPU     = 1, ///< force GPU implementation, try not to store data on CPU when possible
                             HAPI_RES_MGR_DUAL    = 2, ///< create both, duplicate data
};  

struct HAPI_ReserveOpions
{
  int32_t maxGeoms  = 16384;
  int32_t maxImages = 256;
  int32_t maxLights = 131126;
};

HAPI HAPI_SceneLibrary hapiCreateEmpty   (HAPI_RES_MANAGER_TYPE a_type, HAPI_ReserveOpions a_reserveOptions);
HAPI HAPI_SceneLibrary hapiCreateFromFile(HAPI_RES_MANAGER_TYPE a_type, HAPI_ReserveOpions a_reserveOptions, const char* a_filename);
HAPI void              hapiSaveSceneLibrary(HAPI_SceneLibrary, const char* a_filename);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// input types with structures and e.t.c
////

struct HAPI_Geom ///< asbtract type for geometry, usually mesh
{
  int32_t id = -1;
  const char* xmldata = nullptr; ///<! custom fields if we have them
};

struct HAPI_Material ///< asbtract type for material
{
  int32_t id = -1;
  const char* xmldata = nullptr; ///<! custom fields if we have them
};

struct HAPI_Light  ///< asbtract type for light
{
  int32_t id = -1;
  const char* xmldata = nullptr; ///<! custom fields if we have them
};

struct HAPI_Texture ///< asbtract type for images/textures
{
  int32_t id = -1;
  const char* xmldata = nullptr; ///<! custom fields if we have them
};

struct HAPI_Spectrum ///< asbtract type for spectrum
{
  int32_t id = -1;
  const char* xmldata = nullptr; ///<! custom fields if we have them
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum HAPI_CMD_TYPE { HAPI_CMD_APPEND = 0, ///< create new bjects
                     HAPI_CMD_UPDATE = 1, ///< update existing objects
                     HAPI_CMD_DUAL   = 2, ///< simultaniously create and update
                     HAPI_CMD_UNDEFINED = 0xFFFFFFFF,
}; 

struct HAPI_CommandBuffer ///<! use this object to add new data to scene library
{
  int32_t id = -1;
  HAPI_CMD_TYPE type = HAPI_CMD_UNDEFINED;
};

HAPI_CommandBuffer hapiCreateCommandBuffer(HAPI_SceneLibrary a_scnLib, HAPI_CMD_TYPE a_type);
void               hapiCommitCommandBuffer(HAPI_SceneLibrary a_scnLib, HAPI_CommandBuffer a_cmbBuff); ///< Commit and immediately delete it

//// Create new objects 

HAPI_Geom    hapiCreateGeomFromFile   (HAPI_CommandBuffer a_cmdBuff, const char* a_filename);
HAPI_Texture hapiCreateTextureFromFile(HAPI_CommandBuffer a_cmdBuff, const char* a_filename);

//// Update existing objects

HAPI_Geom    hapiGetGeom   (HAPI_CommandBuffer a_cmdBuff, const int32_t a_id);
HAPI_Texture hapiGetTexture(HAPI_CommandBuffer a_cmdBuff, const int32_t a_id);

void         hapiSetGeom   (HAPI_CommandBuffer a_cmdBuff, HAPI_Geom    a_geom);
void         hapiSetTexture(HAPI_CommandBuffer a_cmdBuff, HAPI_Texture a_tex);