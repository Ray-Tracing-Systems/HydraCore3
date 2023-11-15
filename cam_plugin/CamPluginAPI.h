#pragma once
#include <cstdint>      // uint32_t
#include <cstddef>      // for size_t

struct CamParameters    ///<! add any parameter you like to this structure
{
  float fov;
  float aspect;
  float nearPlane;
  float farPlane;
};

static constexpr float CAM_LAMBDA_MIN = 360.0f; ///<! you should statically check that hydra LAMBDA_MIN == CAM_LAMBDA_MIN
static constexpr float CAM_LAMBDA_MAX = 830.0f; ///<! you should statically check that hydra LAMBDA_MAX == CAM_LAMBDA_MAX

struct AuxRayData
{
  uint16_t wavelengthsFixp[4]; ///<! 
  //uint32_t dummy1;           ///<! not used, reserved
  //uint32_t dummy2;           ///<! not used, reserved
};

struct ICamRaysAPI
{
  virtual ~ICamRaysAPI() {}

  /**
   \brief Set camera parameters
   \param a_width         - image width
   \param a_height        - image height
   \param a_camNodeText   - all other camera parameters which you can directly by reading this node
  */
  virtual void SetParameters(int a_width, int a_height, const CamParameters& a_params) {}

  /**
   \brief Put portion of rays in execution queue
   \param out_rayPosAndNear - packed ray origin    (x,y,z) and tNear (w)
   \param out_rayDirAndFar  - packed ray direction (x,y,z) and tFar  (w)
   \param in_blockSize      - ray portion size     (may depend on GPU/device, usually == 1024*512)  
    Please note that it is assumed that rays are uniformly distributed over image plane (and all other integrated dimentions like position on lens)
    for the whole period of time (all passes), the example will be provided.
  */
  virtual void MakeRaysBlock(float* out_rayPosAndNear4f, float* out_rayDirAndFar4f, AuxRayData* out_auxData, size_t in_blockSize, int passId) = 0;

  /**
  \brief Add contribution
  \param out_color  - out float4 image of size a_width*a_height
  \param colors     - in float4 array of size a_size
  \param a_size     - array size
  \param a_width    - image width
  \param a_height   - image height
  */
  virtual void AddSamplesContributionBlock(float* out_color4f, const float* colors4f, size_t in_blockSize, 
                                           uint32_t a_width, uint32_t a_height, int passId) = 0;
};
