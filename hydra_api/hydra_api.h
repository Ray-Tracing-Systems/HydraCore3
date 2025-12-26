// main header for HydraAPI 2.0
#pragma once
#include <cstdint>

#define HR2_SEE_PUGIXML 
#ifdef HR2_SEE_PUGIXML
  #include "pugixml.hpp"
#endif

struct HR2_StorageRef  { int32_t id = -1; }; ///< LEVEL ZERO,  main data storage, scene library, resource manager, navel of the earth. Usually this object is created in a single copy.
struct HR2_SceneRef    { int32_t id = -1; int32_t stgId = -1; }; ///< LEVEL ONE,   per scene data storage. 

struct HR2_GeomRef     { int32_t id = -1; };
struct HR2_MaterialRef { int32_t id = -1; };
struct HR2_LightRef    { int32_t id = -1; };
struct HR2_TextureRef  { int32_t id = -1; };
struct HR2_SpectrumRef { int32_t id = -1; };

struct HR2_CameraRef   { int32_t id = -1; };
struct HR2_SettingsRef { int32_t id = -1; };
struct HR2_FrameImgRef { int32_t id = -1; }; 

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

HR2_StorageRef hr2CreateStorage(HR2_RES_STORAGE_TYPE a_type, HR2_ReserveOpions a_reserveOptions);
HR2_StorageRef hr2CreateStorageFromFile(HR2_RES_STORAGE_TYPE a_type, HR2_ReserveOpions a_reserveOptions, const char* a_filename, bool a_async = false);

void hr2SaveStorage  (HR2_StorageRef a_ref, const char* a_filename, bool a_async = false);
void hr2DeleteStorage(HR2_StorageRef a_ref); ///< detele all
bool hr2StorageIsFinished(HR2_StorageRef a_ref); ///< check whether async scene load/save is completed; use this function within a wait-sleep loop when large scene is loaded/saved
                                                 ///< in the first version async load/save is not planned for implementation, always return true


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum HR2_CMD_TYPE { HR2_CLEAR_AND_APPEND  = 1, ///<! append objects, clear each buffer type that has at least one object in list
                    HR2_UPDATE_AND_APPEND = 2, ///<! update existing objects and append new to the end of each buffer
                    HR2_APPEND_ONLY       = 3, ///<! append objects to the end
                    HR2_UPDATE_ONLY       = 4  ///<! update existing objects only
}; 

enum HR2_CMD_LEVEL { HR2_LVL_STORAGE   = 1, ///<! the most long and heavy operations
                     HR2_LVL_SCENE     = 2  ///<! relatively fast operations
                    }; 

struct HR2_CommandBuffer ///<! use this object to add new data to scene library
{
  int32_t       id    = -1;
  HR2_CMD_TYPE  type  = HR2_CLEAR_AND_APPEND;
  HR2_CMD_LEVEL level = HR2_LVL_STORAGE;
};

HR2_CommandBuffer hr2CommandBufferStorage(HR2_StorageRef a_storage, HR2_CMD_TYPE a_type); ///<! LEVEL ZERO command buffer
HR2_CommandBuffer hr2CommandBufferScene  (HR2_SceneRef   a_scene,   HR2_CMD_TYPE a_type); ///<! LEVEL ONE  command buffer

void              hr2Commit(HR2_CommandBuffer a_cmbBuff, bool a_async = false);           ///<! Commit and then immediately delete it
void              hr2CommitAndRender(HR2_CommandBuffer a_cmbBuff, HR2_CameraRef a_cam, HR2_SettingsRef a_settings, HR2_FrameImgRef a_frameBuffer, bool a_async = false);

/**
\brief these struct is returned by hr2HaveUpdate
*/
struct HR2RenderUpdateInfo
{
  HR2RenderUpdateInfo() : haveUpdateFB(false), haveUpdateMSG(false), finalUpdate(false), progress(0.0f), msg("") {}

  bool        haveUpdateFB;
  bool        haveUpdateMSG;
  bool        finalUpdate;
  float       progress;
  const char* msg;
};

HR2RenderUpdateInfo hr2HaveUpdate(HR2_CommandBuffer a_cmbBuff);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef HR2_SEE_PUGIXML
pugi::xml_node hr2MaterialParamNode(HR2_CommandBuffer a_cmdBuff, HR2_MaterialRef a_mat);
pugi::xml_node hr2LightParamNode   (HR2_CommandBuffer a_cmdBuff, HR2_LightRef    a_lgt);
pugi::xml_node hr2CameraParamNode  (HR2_CommandBuffer a_cmdBuff, HR2_CameraRef   a_cam);
pugi::xml_node hr2SettingsParamNode(HR2_CommandBuffer a_cmdBuff, HR2_SettingsRef a_cam);
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct HR2_FrameBufferInfo
{
  uint32_t width    = 512;
  uint32_t height   = 512;
  uint32_t channels = 4;
  uint32_t bpp      = 16;
};

HR2_MaterialRef hr2CreateMaterial(HR2_CommandBuffer a_cmdBuff);
HR2_LightRef    hr2CreateLight   (HR2_CommandBuffer a_cmdBuff);
HR2_CameraRef   hr2CreateCamera  (HR2_CommandBuffer a_cmdBuff);
HR2_SettingsRef hr2CreateSettings(HR2_CommandBuffer a_cmdBuff);
HR2_SceneRef    hr2CreateScene   (HR2_CommandBuffer a_cmdBuff);
HR2_FrameImgRef hr2CreateFrameImg(HR2_CommandBuffer a_cmdBuff, HR2_FrameBufferInfo a_info);

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
  uint32_t  vertNum    = 0;
  //uint32_t  indicesStride = 3; ///<! 3 for triangle mesh, 4 for quad mesh, quads are not supported yet

  uint32_t* matIdPtr = nullptr;
  uint32_t  matIdAll = 0;
  uint32_t  matIdNum = 1; ///<! if 1, set whole mesh with single material, read matIdAll; else read material indices from matIdPtr
};

HR2_GeomRef    hr2CreateMeshFromData(HR2_CommandBuffer a_cmdBuff, const char* a_meshName, HR2_MeshInput a_input);
HR2_GeomRef    hr2CreateGeomFromFile(HR2_CommandBuffer a_cmdBuff, const char* a_filename);
HR2_TextureRef hr2CreateTextureFromFile(HR2_CommandBuffer a_cmdBuff, const char* a_filename);

int hr2GeomInstance (HR2_CommandBuffer a_cmdBuff, HR2_GeomRef  a_pMesh,  float a_mat[16], const int32_t* a_remapList = nullptr, int32_t a_remapListSize = 0);
int hr2LightInstance(HR2_CommandBuffer a_cmdBuff, HR2_LightRef a_pLight, float a_mat[16], const int32_t* a_remapList = nullptr, int32_t a_remapListSize = 0);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hr2SaveFrameBuffer(HR2_FrameImgRef a_frameImage, const char* a_fileName); 