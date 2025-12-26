#include "hydra_cpu.h"
#include "../integrator_pt.h"

struct CPURenderDriver : HR2::IRenderDriver
{
  CPURenderDriver()
  {
    m_pImpl = std::make_shared<Integrator>(1024*1024, std::vector<uint32_t>());
  }

  void LoadScene(hydra_xml::HydraScene& a_scn) override;
  void CommitDeviceData() override;

  std::shared_ptr<Integrator> m_pImpl;
};

std::shared_ptr<HR2::IRenderDriver> MakeHydraRenderCPU()
{
  return std::make_shared<CPURenderDriver>();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CPURenderDriver::LoadScene(hydra_xml::HydraScene& a_scn)
{
  m_pImpl->LoadScene(a_scn);
}

void CPURenderDriver::CommitDeviceData()
{
  m_pImpl->CommitDeviceData();
}

void HR2::CommandBuffer::CommitToStorage()
{
  // current implementation assume direct forwarding data to renderer
  // m_pRenderImpl = GetRenderDriverForStorage()
  // m_pRenderImpl->LoadScene();
  // m_pRenderImpl->CommitDeviceData();
}