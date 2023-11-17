#pragma once
#include <cstdint>      // uint32_t
#include <cstddef>      // for size_t

struct CamParameters    ///<! add any parameter you like to this structure
{
  float fov;
  float aspect;
  float nearPlane;
  float farPlane;
  int   spectralMode;
};

static constexpr float CAM_LAMBDA_MIN = 360.0f; ///<! you should statically check that hydra LAMBDA_MIN == CAM_LAMBDA_MIN
static constexpr float CAM_LAMBDA_MAX = 830.0f; ///<! you should statically check that hydra LAMBDA_MAX == CAM_LAMBDA_MAX

struct RayPart1 
{
  float    origin[3]; ///<! ray origin, x,y,z
  uint32_t waves01;   ///<! Packed 2 first waves in 16 bit xy, fixed point; 0x0000 => CAM_LAMBDA_MIN; 0xFFFF => CAM_LAMBDA_MAX; wave[0] stored in less significant bits.
};

struct RayPart2 
{
  float    direction[3]; ///<! normalized ray direction, x,y,z
  uint32_t waves23;      ///<! Packed 2 last waves in 16 bit xy, fixed point; 0x0000 => CAM_LAMBDA_MIN; 0xFFFF => CAM_LAMBDA_MAX; wave[1] stored in less significant bits.
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

  virtual void SetBatchSize(int a_tileSize) = 0;

  /**
   \brief Put portion of rays in execution queue
   \param out_rayPosAndNear - packed ray origin    (x,y,z) and tNear (w)
   \param out_rayDirAndFar  - packed ray direction (x,y,z) and tFar  (w)
   \param out_auxData       - packed ray wavelengs and other aux data
   \param in_blockSize      - ray portion size     (may depend on GPU/device, usually == 1024*512)  
    Please note that it is assumed that rays are uniformly distributed over image plane (and all other integrated dimentions like position on lens)
    for the whole period of time (all passes), the example will be provided.
  */
  virtual void MakeRaysBlock(RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f, size_t in_blockSize, int subPassId) = 0;

  /**
  \brief Add contribution
  \param out_color  - out float4 image of size a_width*a_height
  \param colors     - in float4 array of size a_size
  \param a_size     - array size
  \param a_width    - image width
  \param a_height   - image height
  */
  virtual void AddSamplesContributionBlock(float* out_color4f, const float* colors4f, size_t in_blockSize, 
                                           uint32_t a_width, uint32_t a_height, int subPassId) = 0;
};
