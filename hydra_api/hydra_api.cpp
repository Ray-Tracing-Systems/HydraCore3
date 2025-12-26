#include "hydra_api.h"
#include "hydra_cpu.h"
//#include "hydra_vk.h"

#include "LiteScene/hydraxml.h"

#include <memory>
#include <vector>

#include <iostream>
#include <ostream>

struct GlobalContext
{
  static constexpr uint32_t MAX_STORAGES = 4;
  static constexpr uint32_t MAX_COMMMAND_BUFFERS = 8;

  std::shared_ptr<HR2::SceneStorage> storages[MAX_STORAGES] = {};
  uint32_t                           storageTop = 0;

  std::unique_ptr<HR2::CommandBuffer> cmdInFlight[MAX_COMMMAND_BUFFERS];
  std::ostream& textOut = std::cout;
} g_context;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


uint32_t HR2::CommandBuffer::AppendNode(hydra_xml::XML_OBJECT_TYPES a_objType)
{
  const int typeId = int(a_objType);
  assert(typeId < int(hydra_xml::XML_OBJ_TYPES_NUM));
  
  auto nodeNamePair = pStorage->xmlData.RootFor(a_objType);
  assert(nodeNamePair.first != nullptr);

  auto node = nodeNamePair.first.append_child(nodeNamePair.second);
  node.append_attribute(L"id") = pStorage->xmlById[typeId].size();
  pStorage->xmlById[typeId].push_back(node); 
  
  return uint32_t(pStorage->xmlById[typeId].size() - 1);
}

pugi::xml_node HR2::CommandBuffer::NodeById(hydra_xml::XML_OBJECT_TYPES a_objType, uint32_t a_id)
{
  const int typeId = int(a_objType);
  assert(typeId < int(hydra_xml::XML_OBJ_TYPES_NUM));
  return pStorage->xmlById[typeId][a_id];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HR2_StorageRef hr2CreateStorage(HR2_RES_STORAGE_TYPE a_type, HR2_ReserveOpions a_reserveOptions)
{
  HR2_StorageRef res = {};
 
  if(g_context.storageTop >= GlobalContext::MAX_STORAGES)
  {
    g_context.textOut << "[hr2CreateStorage]: MAX_STORAGES exceeded! " << std::endl;
    res.id = -1;
  }
  else
  {
    g_context.storages[g_context.storageTop] = std::make_unique<HR2::SceneStorage>();
    g_context.storages[g_context.storageTop]->xmlData.LoadEmpty();
    res.id = g_context.storageTop;
    g_context.storageTop++;
  }

  return res;
}

HR2_StorageRef hr2CreateStorageFromFile(HR2_RES_STORAGE_TYPE a_type, HR2_ReserveOpions a_reserveOptions, const char* a_filename, bool a_async)
{
  HR2_StorageRef res = {};
  return res;
}

static bool CheckStorage(HR2_StorageRef a_ref, const char* a_funName)
{
  if(a_ref.id >= g_context.storageTop)
  {
    g_context.textOut << "[" << a_funName << "]: Bad storage at: " << a_ref.id << std::endl;
    return false;
  }

  if(g_context.storages[a_ref.id] == nullptr)
  {
    g_context.textOut << "[" << a_funName << "]: Nullptr storage at: " << a_ref.id << std::endl;
    return false;
  }

  return true;
}

static bool CheckCommandBuffer(HR2_CommandBuffer a_ref, const char* a_funName)
{
  if(a_ref.id >= GlobalContext::MAX_COMMMAND_BUFFERS)
  {
    g_context.textOut << "[" << a_funName << "]: Bad Command Buffer at: " << a_ref.id << std::endl;
    return false;
  }

  if(g_context.cmdInFlight[a_ref.id] == nullptr)
  {
    g_context.textOut << "[" << a_funName << "]: Nullptr Command Buffer at: " << a_ref.id << std::endl;
    g_context.cmdInFlight[a_ref.id]->pStorage = nullptr; // always ensure don't have copies
    return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hr2SaveStorage(HR2_StorageRef a_ref, const char* a_filename, bool a_async)
{
  if(!CheckStorage(a_ref, "hr2SaveStorage")) 
    return;
  g_context.storages[a_ref.id]->xmlData.SaveState(a_filename);
}

void hr2DeleteStorage(HR2_StorageRef a_ref)
{
  if(!CheckStorage(a_ref, "hr2DeleteStorage"))
    return; 

  for(uint32_t cmdBuffId=0; cmdBuffId<GlobalContext::MAX_COMMMAND_BUFFERS; cmdBuffId++)
    g_context.cmdInFlight[cmdBuffId] = nullptr;
  g_context.storages[a_ref.id] = nullptr; // unique_ptr should delete it
}

bool hr2StorageIsFinished(HR2_StorageRef a_ref) { return true; }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint32_t FirstEmptyCmdBuff()
{
  uint32_t foundCmdBuffId = uint32_t(-1);
  for(uint32_t cmdBuffId=0; cmdBuffId<GlobalContext::MAX_COMMMAND_BUFFERS; cmdBuffId++)
  {
    if(g_context.cmdInFlight[cmdBuffId] == nullptr)
    {
      foundCmdBuffId = cmdBuffId;
      break;
    }
  }
  return foundCmdBuffId;
}

HR2_CommandBuffer hr2CommandBufferStorage(HR2_StorageRef a_storageRef, HR2_CMD_TYPE a_type)
{
  HR2_CommandBuffer buf = {};
  buf.type  = a_type;
  buf.level = HR2_LVL_STORAGE;
  buf.id    = uint32_t(-1);

  if(!CheckStorage(a_storageRef, "hr2CommandBufferStorage"))
    return buf; 

  const uint32_t foundCmdBuffId = FirstEmptyCmdBuff();
  if(foundCmdBuffId == uint32_t(-1))
  {
    g_context.textOut << "[hr2CommandBufferStorage]: Too many command buffers in flight, max " << GlobalContext::MAX_COMMMAND_BUFFERS << std::endl;
    return buf;
  }
  
  buf.id = foundCmdBuffId;
  g_context.cmdInFlight[buf.id] = std::make_unique<HR2::CommandBuffer>(g_context.storages[a_storageRef.id]);
  g_context.cmdInFlight[buf.id]->m_type  = buf.type;
  g_context.cmdInFlight[buf.id]->m_level = buf.level;
  g_context.cmdInFlight[buf.id]->m_stgId = a_storageRef.id;

  return buf;
}

HR2_CommandBuffer hr2CommandBufferScene(HR2_SceneRef a_scene,  HR2_CMD_TYPE a_type)
{
  HR2_CommandBuffer buf = {};
  buf.type  = a_type;
  buf.level = HR2_LVL_SCENE;

  if(a_scene.stgId >= g_context.MAX_STORAGES)
  {
    g_context.textOut << "[hr2CommandBufferScene]: Bad scene storage id " << a_scene.stgId << ", scene: " << a_scene.id << std::endl;
    return buf;
  }

  //TODO: validate scene was actually created!

  const uint32_t foundCmdBuffId = FirstEmptyCmdBuff();
  if(foundCmdBuffId == uint32_t(-1))
  {
    g_context.textOut << "[hr2CommandBufferScene]: Too many command buffers in flight, max " << GlobalContext::MAX_COMMMAND_BUFFERS << std::endl;
    return buf;
  }

  buf.id = foundCmdBuffId;
  g_context.cmdInFlight[buf.id] = std::make_unique<HR2::CommandBuffer>(g_context.storages[a_scene.stgId]);
  g_context.cmdInFlight[buf.id]->m_type  = buf.type;
  g_context.cmdInFlight[buf.id]->m_level = buf.level;
  g_context.cmdInFlight[buf.id]->m_stgId = a_scene.stgId;

  return buf;
}

void hr2Commit(HR2_CommandBuffer a_cmbBuff, bool a_async)
{
  if(!CheckCommandBuffer(a_cmbBuff, "hr2Commit"))
    return;
  
  g_context.cmdInFlight[a_cmbBuff.id]->CommitToStorage();  
  g_context.cmdInFlight[a_cmbBuff.id] = nullptr;
}

void hr2CommitAndRender(HR2_CommandBuffer a_cmbBuff, HR2_CameraRef a_cam, HR2_SettingsRef a_settings, HR2_FrameImgRef a_frameBuffer, bool a_async)
{
  hr2Commit(a_cmbBuff, a_async);
  //todo: call actual render
}

HR2RenderUpdateInfo hr2HaveUpdate(HR2_CommandBuffer a_cmbBuff)
{
  HR2RenderUpdateInfo res{};
  return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HR2_GeomRef hr2CreateMeshFromData(HR2_CommandBuffer a_cmdBuff, const char* a_meshName, HR2_MeshInput a_input)
{
  HR2_GeomRef res = {};
  res.id = g_context.cmdInFlight[a_cmdBuff.id]->AppendNode(hydra_xml::XML_OBJ_GEOMETRY);
           g_context.cmdInFlight[a_cmdBuff.id]->meshPtrById[res.id] = a_input;
  
  auto node = g_context.cmdInFlight[a_cmdBuff.id]->NodeById(hydra_xml::XML_OBJ_GEOMETRY, res.id);
  node.append_attribute(L"ptrs") = 1;
  return res;
}

HR2_MaterialRef hr2CreateMaterial(HR2_CommandBuffer a_cmdBuff)
{
  HR2_MaterialRef res = {};
  res.id = g_context.cmdInFlight[a_cmdBuff.id]->AppendNode(hydra_xml::XML_OBJ_MATERIALS);
  return res;
}

HR2_LightRef  hr2CreateLight(HR2_CommandBuffer a_cmdBuff)
{
  HR2_LightRef res = {};
  res.id = g_context.cmdInFlight[a_cmdBuff.id]->AppendNode(hydra_xml::XML_OBJ_LIGHT);
  return res;
}

HR2_CameraRef  hr2CreateCamera(HR2_CommandBuffer a_cmdBuff)
{
  HR2_CameraRef res = {};
  res.id = g_context.cmdInFlight[a_cmdBuff.id]->AppendNode(hydra_xml::XML_OBJ_CAMERA);
  return res;
}

HR2_SettingsRef hr2CreateSettings(HR2_CommandBuffer a_cmdBuff)
{
  HR2_SettingsRef res = {};
  res.id = g_context.cmdInFlight[a_cmdBuff.id]->AppendNode(hydra_xml::XML_OBJ_SETTINGS);
  return res;
}

HR2_SceneRef    hr2CreateScene   (HR2_CommandBuffer a_cmdBuff)
{
  HR2_SceneRef res = {};
  res.id    = g_context.cmdInFlight[a_cmdBuff.id]->AppendNode(hydra_xml::XML_OBJ_SCENE);
  res.stgId = g_context.cmdInFlight[a_cmdBuff.id]->m_stgId;
  return res;
}

HR2_FrameImgRef hr2CreateFrameImg(HR2_CommandBuffer a_cmdBuff, HR2_FrameBufferInfo a_info)
{
  HR2_FrameImgRef res = {};
  return res;
}

#ifdef HR2_SEE_PUGIXML

pugi::xml_node hr2MaterialParamNode(HR2_CommandBuffer a_cmdBuff, HR2_MaterialRef a_mat)
{
  return g_context.cmdInFlight[a_cmdBuff.id]->NodeById(hydra_xml::XML_OBJ_MATERIALS, a_mat.id);
}

pugi::xml_node hr2LightParamNode(HR2_CommandBuffer a_cmdBuff, HR2_LightRef a_light)
{
  return g_context.cmdInFlight[a_cmdBuff.id]->NodeById(hydra_xml::XML_OBJ_LIGHT, a_light.id);
}

pugi::xml_node hr2CameraParamNode(HR2_CommandBuffer a_cmdBuff, HR2_CameraRef a_cam)
{
  return g_context.cmdInFlight[a_cmdBuff.id]->NodeById(hydra_xml::XML_OBJ_CAMERA, a_cam.id);
}

pugi::xml_node hr2SettingsParamNode(HR2_CommandBuffer a_cmdBuff, HR2_SettingsRef a_settings)
{
  return g_context.cmdInFlight[a_cmdBuff.id]->NodeById(hydra_xml::XML_OBJ_SETTINGS, a_settings.id);
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int hr2GeomInstance (HR2_CommandBuffer a_cmdBuff, HR2_GeomRef  a_pMesh,  float a_mat[16], const int32_t* a_remapList, int32_t a_remapListSize)
{
  return 0;
}

int hr2LightInstance(HR2_CommandBuffer a_cmdBuff, HR2_LightRef a_pLight, float a_mat[16], const int32_t* a_remapList, int32_t a_remapListSize)
{
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hr2SaveFrameBuffer(HR2_FrameImgRef a_frameImage, const char* a_fileName)
{

}
