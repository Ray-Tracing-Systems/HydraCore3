#include "hydra_api.h"
#include "hydra_cpu.h"
//#include "hydra_vk.h"

//struct Impl_XXX
//{
//  CPU_IMPL* pImplCPU;
//  CPU_IMPL* pImplGPU;
//};

HAPI HR2_SceneLibraryRef hr2CreateLibrary(HR2_RES_STORAGE_TYPE a_type, HR2_ReserveOpions a_reserveOptions)
{
  HR2_SceneLibraryRef res = {};
  return res;
}

HAPI HR2_CommandBuffer hr2CreateCommandBuffer(HR2_SceneLibraryRef a_scnLib, HR2_CMD_TYPE a_type)
{
  HR2_CommandBuffer buf = {};
  return buf;
}

HAPI void hr2CommitCommandBuffer(HR2_CommandBuffer a_cmbBuff, bool a_async)
{

}

HAPI HR2_GeomRef hr2CreateMeshFromData(HR2_CommandBuffer a_cmdBuff, const char* a_meshName, HR2_MeshInput a_input)
{
  HR2_GeomRef res = {};
  return res;
}

HAPI HR2_MaterialRef hr2CreateMaterial(HR2_CommandBuffer a_cmdBuff, const char* a_matName)
{
  HR2_MaterialRef res = {};
  return res;
}

HAPI HR2_LightRef  hr2CreateLight(HR2_CommandBuffer a_cmdBuff, const char* a_lgtName)
{
  HR2_LightRef res = {};
  return res;
}

HAPI HR2_CameraRef  hr2CreateCamera(HR2_CommandBuffer a_cmdBuff, const char* a_camName)
{
  HR2_CameraRef res = {};
  return res;
}

#ifdef HAPI_SEE_PUGIXML

HAPI pugi::xml_node hr2MaterialParamNode(HR2_MaterialRef a_mat)
{
  return pugi::xml_node();
}

HAPI pugi::xml_node hr2LightParamNode(HR2_LightRef a_mat)
{
  return pugi::xml_node();
}

HAPI pugi::xml_node hr2CameraParamNode(HR2_CameraRef a_cam)
{
  return pugi::xml_node();
}

#endif