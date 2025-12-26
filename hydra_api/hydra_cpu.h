// CPU implementation of HydraAPI 2.0 resource manager
#pragma once

#include "hydra_api.h"
#include "LiteScene/hydraxml.h"

#include <memory>

namespace HR2
{
  struct IRenderDriver
  {
    virtual void LoadScene(hydra_xml::HydraScene& a_scn) = 0;
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
    

    uint32_t       AppendNode(hydra_xml::XML_OBJECT_TYPES a_objType);
    pugi::xml_node NodeById  (hydra_xml::XML_OBJECT_TYPES a_objType, uint32_t a_id);

    virtual void CommitToStorage();
  };

};

