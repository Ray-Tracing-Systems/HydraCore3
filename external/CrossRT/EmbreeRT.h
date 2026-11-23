#pragma once

#include <vector>
#include <cstdint>

#include "CrossRT.h"
#include "embree3/rtcore.h"

/**
\brief CPU implementation of ISceneObject via Embree; obtain it with dynamic_cast<EmbreeRT*>(pScene->UnderlyingImpl(0))
*/
class EmbreeRT : public ISceneObject
{
public:
  EmbreeRT();
  ~EmbreeRT();
  void ClearGeom() override;
  const char* Name() const override { return "VulkanRTX"; }

  uint32_t AddGeom_Triangles3f(const float * a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, 
                               uint32_t a_flags, size_t vByteStride) override;
  void     UpdateGeom_Triangles3f(uint32_t a_geomId, const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices,
                                  size_t a_indNumber, uint32_t a_flags, size_t vByteStride) override;

  void ClearScene() override; 
  void CommitScene  (uint32_t options = BUILD_HIGH) override; 
  
  uint32_t AddInstance(uint32_t a_geomId, const LiteMath::float4x4& a_matrix) override;
  uint32_t AddInstanceMotion(uint32_t a_geomId, const LiteMath::float4x4* a_matrices, uint32_t a_matrixNumber) override;
  void     UpdateInstance(uint32_t a_instanceId, const LiteMath::float4x4& a_matrix) override;

  CRT_Hit  RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) override;
  bool     RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) override;
  CRT_Hit  RayQuery_NearestHitMotion(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, float time) override;
  bool     RayQuery_AnyHitMotion(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, float time) override;

  // direct access to Embree data for code which calls Embree itself (ISPC back end of kernel_slicer);
  // CRT_Hit::geomId of a hit is GetGeomIdByInstId()[rtcHit.instID[0]]
  //
  RTCScene        GetEmbreeScene()    const { return m_scene; }
  const uint32_t* GetGeomIdByInstId() const { return m_geomIdByInstId.data(); }

protected:
  uint32_t maxFrames = 1;
  RTCDevice m_device = nullptr;
  RTCScene  m_scene  = nullptr;

  struct BlasInfo
  {
    size_t verticesNum = 0;
    size_t indicesNum  = 0;
    size_t vByteStride = 0;
  };

  std::vector<RTCScene>    m_blas;
  std::vector<BlasInfo>    m_blasInfo;
  std::vector<RTCGeometry> m_inst;
  std::vector<uint32_t>    m_geomIdByInstId;
};
