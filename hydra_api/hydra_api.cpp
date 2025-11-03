#include "hydra_api.h"
#include "hydra_cpu.h"
//#include "hydra_vk.h"

//struct Impl_XXX
//{
//  CPU_IMPL* pImplCPU;
//  CPU_IMPL* pImplGPU;
//};

HAPI HAPI_SceneLibrary hapiCreateLibraryEmpty(HAPI_RES_STORAGE_TYPE a_type, HAPI_ReserveOpions a_reserveOptions)
{
  HAPI_SceneLibrary res = {};
  return res;
}

HAPI HAPI_CommandBuffer hapiCreateCommandBuffer(HAPI_SceneLibrary a_scnLib, HAPI_CMD_TYPE a_type)
{
  HAPI_CommandBuffer buf = {};
  return buf;
}

HAPI void hapiCommitCommandBuffer(HAPI_CommandBuffer a_cmbBuff, bool a_async)
{

}

HAPI HAPI_Geom hapiCreateMeshFromData(HAPI_CommandBuffer a_cmdBuff, const char* a_meshName, HAPI_MeshInput a_input)
{
  HAPI_Geom res = {};
  return res;
}

HAPI HAPI_Material hapiCreateMaterial(HAPI_CommandBuffer a_cmdBuff, const char* a_matName)
{
  HAPI_Material res = {};
  return res;
}

HAPI HAPI_Light  hapiCreateLight(HAPI_CommandBuffer a_cmdBuff, const char* a_lgtName)
{
  HAPI_Light res = {};
  return res;
}

HAPI HAPI_Camera  hapiCreateCamera(HAPI_CommandBuffer a_cmdBuff, const char* a_camName)
{
  HAPI_Camera res = {};
  return res;
}

#ifdef HAPI_SEE_PUGIXML

HAPI pugi::xml_node hapiMaterialParamNode(HAPI_Material a_mat)
{
  return pugi::xml_node();
}

HAPI pugi::xml_node hapiLightParamNode(HAPI_Light a_mat)
{
  return pugi::xml_node();
}

HAPI pugi::xml_node hapiCameraParamNode(HAPI_Camera a_cam)
{
  return pugi::xml_node();
}

#endif