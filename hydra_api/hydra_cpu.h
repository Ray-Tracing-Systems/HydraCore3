// CPU implementation of HydraAPI 2.0 resource manager
#pragma once

#include "hydra_api.h"
#include "LiteScene/hydraxml.h"

static constexpr uint32_t SCN_UPDATE_TEXTURES    = 1;
static constexpr uint32_t SCN_UPDATE_SPECTRUM    = 2;
static constexpr uint32_t SCN_UPDATE_LIGHTS      = 4;
static constexpr uint32_t SCN_UPDATE_MATERIALS   = 8;
static constexpr uint32_t SCN_UPDATE_CAMERA      = 16;
static constexpr uint32_t SCN_UPDATE_GEOMETRY    = 32;
static constexpr uint32_t SCN_UPDATE_INSTANCES   = 64;
static constexpr uint32_t SCN_UPDATE_REMAP_LISTS = 128;
static constexpr uint32_t SCN_UPDATE_SETTINGS    = 256;
static constexpr uint32_t SCN_UPDATE_ALL         = 0xFFFFFFFF;

#include <memory>

namespace HR2
{
  struct RDScene_Input
  {
    const std::unordered_map<int, HR2_MeshInput>* pMeshPtrs = nullptr;
  };

  struct IRenderDriver
  {
    IRenderDriver(){}
    virtual ~IRenderDriver(){}
    virtual void LoadScene(hydra_xml::HydraScene& a_scn, const RDScene_Input& a_input, uint32_t a_updateFlags) = 0;
    virtual void CommitDeviceData() = 0;
  };
  
  std::shared_ptr<IRenderDriver> MakeHydraRenderCPU(); 

  struct SceneStorage
  {
    SceneStorage() {
      for(int i=0;i<hydra_xml::XML_OBJ_TYPES_NUM; i++)
        xmlById[i].reserve(1024); 

      m_pDriver = MakeHydraRenderCPU();
    }
    virtual ~SceneStorage(){}

    hydra_xml::HydraScene       xmlData;
    std::vector<pugi::xml_node> xmlById[hydra_xml::XML_OBJ_TYPES_NUM]; 

    std::shared_ptr<IRenderDriver> m_pDriver;
  };

  struct CommandBuffer
  {
    CommandBuffer(){}
    CommandBuffer(std::shared_ptr<SceneStorage> a_pStorage) : pStorage(a_pStorage) {}
    virtual ~CommandBuffer(){}

    std::shared_ptr<SceneStorage> pStorage = nullptr;
    int32_t                       m_stgId  = -1;
    HR2_CMD_TYPE                  m_type   = HR2_APPEND_ONLY;
    HR2_CMD_LEVEL                 m_level  = HR2_LVL_SCENE;
    uint32_t                      m_updateFlags = 0;

    uint32_t       AppendNode(hydra_xml::XML_OBJECT_TYPES a_objType);
    pugi::xml_node NodeById  (hydra_xml::XML_OBJECT_TYPES a_objType, uint32_t a_id);

    virtual void CommitToStorage();

    /// @brief /// temporary storage for mesh pointers by mesh id to pass them inside render
    std::unordered_map<int, HR2_MeshInput> meshPtrById;
    
  };

};

