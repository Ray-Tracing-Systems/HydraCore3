// CPU implementation of HydraAPI 2.0 resource manager
#pragma once

#include "hydra_api.h"
#include "LiteScene/hydraxml.h"

#include <memory>

namespace HR2
{
  static constexpr uint32_t MAX_FRAME_IMAGES = 8;

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

    virtual void Render(uint32_t startX, uint32_t startY, int32_t sizeX, uint32_t sizeY, uint32_t channels, float* data, uint32_t a_passNumber) = 0;
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
    
    std::vector<float> fbData[MAX_FRAME_IMAGES];
    uint3              fbSize[MAX_FRAME_IMAGES];
    uint32_t           fbTop = 0;
  };

  struct CommandBuffer
  {
    CommandBuffer(){}
    CommandBuffer(std::shared_ptr<SceneStorage> a_pStorage) : pStorage(a_pStorage) {}
    virtual ~CommandBuffer(){}

    std::shared_ptr<SceneStorage> pStorage = nullptr;
    int32_t                       m_stgId  = -1;
    int32_t                       m_scnId  = -1;
    HR2_CMD_TYPE                  m_type   = HR2_APPEND_ONLY;
    HR2_CMD_LEVEL                 m_level  = HR2_LVL_SCENE;
    uint32_t                      m_updateFlags = 0;

    uint32_t       AppendNode(hydra_xml::XML_OBJECT_TYPES a_objType);
    pugi::xml_node NodeById  (hydra_xml::XML_OBJECT_TYPES a_objType, uint32_t a_id);
    
    virtual void CommitToStorage();

    /// @brief /// temporary storage for mesh pointers by mesh id to pass them inside render
    std::unordered_map<int, HR2_MeshInput> meshPtrById;

    /// @brief /// scene node, instances and lights
    pugi::xml_node m_sceneNode;
    int32_t instTop = 0;
    int32_t lghtTop = 0;
  };

};

