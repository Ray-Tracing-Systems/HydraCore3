#include "hydra_cpu.h"
#include "../integrator_pt.h"

struct HydraCore3RenderDriver : HR2::IRenderDriver
{
  HydraCore3RenderDriver()
  {
    m_pImpl = std::make_shared<Integrator>(1024*1024, std::vector<uint32_t>());
  }

  void LoadScene(hydra_xml::HydraScene& a_scn) override;
  void CommitDeviceData() override;

  std::shared_ptr<Integrator> m_pImpl;
};

std::shared_ptr<HR2::IRenderDriver> HR2::MakeHydraRenderCPU()
{
  return std::make_shared<HydraCore3RenderDriver>();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HydraCore3RenderDriver::LoadScene(hydra_xml::HydraScene& a_scn)
{
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

  pStorage->m_pDriver->LoadScene(pStorage->xmlData);
  pStorage->m_pDriver->CommitDeviceData();
}