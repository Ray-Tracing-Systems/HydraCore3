#include "external/HeavyRT/core/CrossRT.h"
#include <algorithm>

// HeavyRT's FactoryRT.cpp also pulls in every renderer implementation. Hydra
// instantiates BVH2FatRT directly, and only needs these common metric methods.
MetricStats ISceneObject::GetStats()
{
  MetricStats stats = {};
  const double rays = std::max(double(m_stats.raysNumber), 1.0);
  stats.avgNC = float(double(m_stats.NC) / rays);
  stats.avgLC = float(double(m_stats.LC) / rays);
  stats.avgTC = float(double(m_stats.TC) / rays);
  for(int i = 0; i < TREELET_ARR_SIZE; ++i)
  {
    stats.avgLJC[i] = float(double(m_stats.LJC[i]) / rays);
    stats.avgCMC[i] = float(double(m_stats.CMC[i]) / rays);
    stats.avgWSS[i] = uint32_t(m_stats.WSS[i].size());
  }
  stats.avgBLB = float(double(m_stats.BLB) / rays);
  stats.avgSOC = float(double(m_stats.SOC) / rays);
  stats.avgSBL = float(double(m_stats.SBL) / rays);
  stats.bvhTotalSize = m_stats.bvhTotalSize;
  stats.geomTotalSize = m_stats.geomTotalSize;
  return stats;
}

void ISceneObject::ResetStats() { m_stats.clear(); }
