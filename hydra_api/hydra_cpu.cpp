#include "hydra_cpu.h"
#include "../integrator_pt.h"

struct HydraCore3RenderDriver : HR2::IRenderDriver
{
  HydraCore3RenderDriver()
  {
    m_pImpl = std::make_shared<Integrator>(1024*1024, std::vector<uint32_t>());
  }

  void LoadScene(hydra_xml::HydraScene& a_scn, const HR2::RDScene_Input& a_input) override;
  void CommitDeviceData() override;

  std::shared_ptr<Integrator> m_pImpl;
};

std::shared_ptr<HR2::IRenderDriver> HR2::MakeHydraRenderCPU()
{
  return std::make_shared<HydraCore3RenderDriver>();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HydraCore3RenderDriver::LoadScene(hydra_xml::HydraScene& a_scn, const HR2::RDScene_Input& a_input)
{
  // put / convert mesh  pointers / data, then pass data to LoadScene
  //
  std::unordered_map<int, Mesh4fInput> meshPtrs(a_input.pMeshPtrs->size());
  
  for(auto p = a_input.pMeshPtrs->begin(); p != a_input.pMeshPtrs->end(); ++p)
  {
    Mesh4fInput input;
    input.vPosPtr        = p->second.vPosPtr;
    input.vPosByteStride = p->second.vPosStride*sizeof(float);
    input.vNormPtr4f     = p->second.vNormPtr;
    input.vTangPtr4f     = p->second.vTangPtr;
    input.vTexCoord2f    = p->second.vTexCoordPtr;
    input.vertNum        = p->second.vertNum;
    
    input.indicesPtr = p->second.indicesPtr;
    input.indicesNum = p->second.indicesNum;
    
    input.matIdPtr   = p->second.matIdPtr;
    input.matIdAll   = p->second.matIdAll;
    input.matIdNum   = p->second.matIdNum;
    
    if(p->second.vNormStride != 4)
    {
      std::cout << "[HydraCore3RenderDriver::LoadScene]: unsuported normal stride = " << p->second.vNormStride << std::endl;
      exit(0);
    }

    if(p->second.vTangStride != 4)
    {
      std::cout << "[HydraCore3RenderDriver::LoadScene]: unsuported tangent stride = " << p->second.vTangStride  << std::endl;
      exit(0);
    }

    if(p->second.vTexCoordStride != 2)
    {
      std::cout << "[HydraCore3RenderDriver::LoadScene]: unsuported texcoord stride = " << p->second.vTexCoordPtr << std::endl;
      exit(0);
    }
    
    meshPtrs[p->first] = input;
  }

  m_pImpl->LoadScene_SetMeshPointers(&meshPtrs);

  // put / convert image pointers / data, then pass data to LoadScene
  //

  //m_pImpl->LoadScene_SetImagePointers(...);

  m_pImpl->LoadScene(a_scn);
}

void HydraCore3RenderDriver::CommitDeviceData()
{
  m_pImpl->CommitDeviceData();
}

void HR2::CommandBuffer::CommitToStorage()
{
  if(pStorage == nullptr)
    return;
  if(pStorage->m_pDriver == nullptr)
    return;

  pStorage->xmlData.SaveState("z_debug.xml");
  
  HR2::RDScene_Input input;
  input.pMeshPtrs = &meshPtrById;
  //input.pImagePtrs = 

  pStorage->m_pDriver->LoadScene(pStorage->xmlData, input);
  pStorage->m_pDriver->CommitDeviceData();
}