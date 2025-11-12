#include "hydra_api.h"
#include "hydra_cpu.h"
//#include "hydra_vk.h"

//struct Impl_XXX
//{
//  CPU_IMPL* pImplCPU;
//  CPU_IMPL* pImplGPU;
//};

HR2_SceneLibraryRef hr2CreateLibrary(HR2_RES_STORAGE_TYPE a_type, HR2_ReserveOpions a_reserveOptions)
{
  HR2_SceneLibraryRef res = {};
  return res;
}

HR2_CommandBuffer hr2CreateCommandBuffer(HR2_SceneLibraryRef a_scnLib, HR2_CMD_TYPE a_type, HR2_CMD_LEVEL a_lvl)
{
  HR2_CommandBuffer buf = {};
  buf.type  = a_type;
  buf.level = a_lvl;
  return buf;
}

void hr2CommitCommandBuffer(HR2_CommandBuffer a_cmbBuff, bool a_async)
{

}

HR2_GeomRef hr2CreateMeshFromData(HR2_CommandBuffer a_cmdBuff, const char* a_meshName, HR2_MeshInput a_input)
{
  HR2_GeomRef res = {};
  return res;
}

HR2_MaterialRef hr2CreateMaterial(HR2_CommandBuffer a_cmdBuff)
{
  HR2_MaterialRef res = {};
  return res;
}

HR2_LightRef  hr2CreateLight(HR2_CommandBuffer a_cmdBuff)
{
  HR2_LightRef res = {};
  return res;
}

HR2_CameraRef  hr2CreateCamera(HR2_CommandBuffer a_cmdBuff)
{
  HR2_CameraRef res = {};
  return res;
}

HR2_SettingsRef hr2CreateSettings(HR2_CommandBuffer a_cmdBuff)
{
  HR2_SettingsRef res = {};
  return res;
}

HR2_SceneRef    hr2CreateScene   (HR2_CommandBuffer a_cmdBuff)
{
  HR2_SceneRef res = {};
  return res;
}

HR2_FrameImgRef hr2CreateFrameImg(HR2_CommandBuffer a_cmdBuff, HR2_FrameBufferInfo a_info)
{
  HR2_FrameImgRef res = {};
  return res;
}

#ifdef HR2_SEE_PUGIXML

pugi::xml_node hr2MaterialParamNode(HR2_MaterialRef a_mat)
{
  return pugi::xml_node();
}

pugi::xml_node hr2LightParamNode(HR2_LightRef a_mat)
{
  return pugi::xml_node();
}

pugi::xml_node hr2CameraParamNode(HR2_CameraRef a_cam)
{
  return pugi::xml_node();
}

pugi::xml_node hr2SettingsParamNode(HR2_SettingsRef a_cam)
{
  return pugi::xml_node();
}

#endif