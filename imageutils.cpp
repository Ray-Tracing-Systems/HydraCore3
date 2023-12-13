#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#if defined(_WIN32)
    #ifndef NOMINMAX
    #define NOMINMAX
    #endif
#endif

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

bool SaveImage4fToEXR(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, bool a_invertY = false) 
{
  EXRHeader header;
  InitEXRHeader(&header);

  EXRImage image;
  InitEXRImage(&image);
  image.num_channels = 3;

  std::vector<float> images[3];
  images[0].resize(width * height);
  images[1].resize(width * height);
  images[2].resize(width * height);

  // Split RGBARGBARGBA... into R, G and B layer
  if(a_invertY) {
    for(int y=0;y<height;y++) {
      const int offsetY1 = y*width*4;
      const int offsetY2 = (height-y-1)*width*4;
      for(int x=0;x<width;x++) {
        images[0][(offsetY1 >> 2) + x] = rgb[offsetY2 + x*4 + 0]*a_normConst;
        images[1][(offsetY1 >> 2) + x] = rgb[offsetY2 + x*4 + 1]*a_normConst;
        images[2][(offsetY1 >> 2) + x] = rgb[offsetY2 + x*4 + 2]*a_normConst; 
      }
    }   
  }
  else {
    for (size_t i = 0; i < size_t(width * height); i++) {
      images[0][i] = rgb[4*i+0]*a_normConst;
      images[1][i] = rgb[4*i+1]*a_normConst;
      images[2][i] = rgb[4*i+2]*a_normConst;
    }
  }

  float* image_ptr[3];
  image_ptr[0] = images[2].data(); // B
  image_ptr[1] = images[1].data(); // G
  image_ptr[2] = images[0].data(); // R

  image.images = (unsigned char**)image_ptr;
  image.width  = width;
  image.height = height;
  header.num_channels = 3;
  header.channels     = (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels);
  // Must be (A)BGR order, since most of EXR viewers expect this channel order.
  strncpy(header.channels[0].name, "B", 255); header.channels[0].name[strlen("B")] = '\0';
  strncpy(header.channels[1].name, "G", 255); header.channels[1].name[strlen("G")] = '\0';
  strncpy(header.channels[2].name, "R", 255); header.channels[2].name[strlen("R")] = '\0';

  header.pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
  header.requested_pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
  for (int i = 0; i < header.num_channels; i++) {
    header.pixel_types[i]           = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
    header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of output image to be stored in .EXR
  }
 
  const char* err = nullptr; 
  int ret = SaveEXRImageToFile(&image, &header, outfilename, &err);
  if (ret != TINYEXR_SUCCESS) {
    fprintf(stderr, "Save EXR err: %s\n", err);
    FreeEXRErrorMessage(err); // free's buffer for an error message
    return false;
  }
  //printf("Saved exr file. [%s] \n", outfilename);

  free(header.channels);
  free(header.pixel_types);
  free(header.requested_pixel_types);

  return true;
}

void FlipYAndNormalizeImage3D(float* data, int width, int height, int channels, float a_normConst = 1.0f)
{
  const size_t size = size_t(width*height*channels);
  for(size_t i=0;i<size;i++)
    data[i] *= a_normConst;
  
  std::vector<float> line(width*channels);
  for(int y1=0;y1<height;y1++) 
  {
    const int y2 = height-y1-1;  
    if(y1 == y2)
      continue;
  
    const int offsetY1 = y1*width*channels;
    const int offsetY2 = y2*width*channels; 
    
    memcpy(line.data(),     data + offsetY1, width*channels*sizeof(float));
    memcpy(data + offsetY1, data + offsetY2, width*channels*sizeof(float));
    memcpy(data + offsetY2, line.data(),     width*channels*sizeof(float));
  }
}

bool SaveImage3DToEXR(const float* data, int width, int height, int channels, const char* outfilename) 
{
  EXRHeader header;
  InitEXRHeader(&header);

  EXRImage image;
  InitEXRImage(&image);
  image.num_channels = 1;

  const float* image_ptr[1] = {data};

  image.images = (unsigned char**)image_ptr;
  image.width  = width;
  image.height = height;
  header.num_channels = 1;
  header.channels     = (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels);

  strncpy(header.channels[0].name, "Y", 255); header.channels[0].name[strlen("Y")] = '\0';

  header.pixel_types           = (int*)malloc(sizeof(int) * header.num_channels);
  header.requested_pixel_types = (int*)malloc(sizeof(int) * header.num_channels);
  for (int i = 0; i < header.num_channels; i++) {
    header.pixel_types[i]           = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
    header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of output image to be stored in .EXR
  }
 
  const char* err = nullptr; 
  int ret = SaveEXRImageToFile(&image, &header, outfilename, &err);
  if (ret != TINYEXR_SUCCESS) {
    fprintf(stderr, "Save EXR err: %s\n", err);
    FreeEXRErrorMessage(err); // free's buffer for an error message
    return false;
  }
  //printf("Saved exr file. [%s] \n", outfilename);

  free(header.channels);
  free(header.pixel_types);
  free(header.requested_pixel_types);

  return true;
}

void FlipYAndSaveFrameBufferToEXR(float* data, int width, int height, int channels, const char* outfilename, float a_normConst = 1.0f)
{
  if(channels == 4)
    SaveImage4fToEXR(data, width, height, outfilename, a_normConst, true);
  else
  {
    //FlipYAndNormalizeImage3D(data, width, height, channels, a_normConst);
    std::vector<float> image2(width*height*channels);
    for(int y=0;y<height;y++) {
      const int offsetY1 = y*width*channels;
      const int offsetY2 = (height-y-1)*width*channels;
      for(int x=0;x<width*channels;x++) {
        image2[offsetY2+x] = data[offsetY1+x]*a_normConst;
      }
    }   
    SaveImage3DToEXR(image2.data(), width, height, channels, outfilename);
  }
}

static inline float clamp(float u, float a, float b) { return std::min(std::max(a, u), b); }

static inline unsigned RealColorToUint32(float real_color[4])
{
  float  r = real_color[0]*255.0f;
  float  g = real_color[1]*255.0f;
  float  b = real_color[2]*255.0f;
  float  a = real_color[3]*255.0f;

  unsigned char red   = (unsigned char)r;
  unsigned char green = (unsigned char)g;
  unsigned char blue  = (unsigned char)b;
  unsigned char alpha = (unsigned char)a;

  return red | (green << 8) | (blue << 16) | (alpha << 24);
}

struct Pixel { unsigned char r, g, b; };

static bool WriteBMP(const char* fname, Pixel* a_pixelData, int width, int height)
{
  int paddedsize = (width*height) * sizeof(Pixel);
  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};

  bmpfileheader[ 2] = (unsigned char)(paddedsize    );
  bmpfileheader[ 3] = (unsigned char)(paddedsize>> 8);
  bmpfileheader[ 4] = (unsigned char)(paddedsize>>16);
  bmpfileheader[ 5] = (unsigned char)(paddedsize>>24);

  bmpinfoheader[ 4] = (unsigned char)(width    );
  bmpinfoheader[ 5] = (unsigned char)(width>> 8);
  bmpinfoheader[ 6] = (unsigned char)(width>>16);
  bmpinfoheader[ 7] = (unsigned char)(width>>24);
  bmpinfoheader[ 8] = (unsigned char)(height    );
  bmpinfoheader[ 9] = (unsigned char)(height>> 8);
  bmpinfoheader[10] = (unsigned char)(height>>16);
  bmpinfoheader[11] = (unsigned char)(height>>24);

  std::ofstream out(fname, std::ios::out | std::ios::binary);
  if(!out.is_open())
    return false;
  out.write((const char*)bmpfileheader, 14);
  out.write((const char*)bmpinfoheader, 40);
  out.write((const char*)a_pixelData, paddedsize);
  out.flush();
  out.close();
  return true;
}

static bool SaveBMP(const char* fname, const unsigned int* pixels, int w, int h)
{
  std::vector<Pixel> pixels2(w*h);

  for (size_t i = 0; i < pixels2.size(); i++)
  {
    Pixel px;
    px.r       = (pixels[i] & 0x00FF0000) >> 16;
    px.g       = (pixels[i] & 0x0000FF00) >> 8;
    px.b       = (pixels[i] & 0x000000FF);
    pixels2[i] = px;
  }

  return WriteBMP(fname, &pixels2[0], w, h);
}

static inline float linearToSRGB(float l)
{
  if(l <= 0.00313066844250063f)
    return l * 12.92f;
  else
    return 1.055f * std::pow(l, 1.0f/2.4f) - 0.055f;
}

std::vector<uint32_t> FrameBufferColorToLDRImage(const float* rgb, int width, int height, float a_normConst, float a_gamma)
{
  std::vector<uint32_t> pixelData(width*height);
  if(std::abs(a_gamma - 2.4f) < 1e-5f) {
    for(int i=0;i<width*height;i++)
    {
      float color[4];
      color[0]     = linearToSRGB(clamp(rgb[4*i+0]*a_normConst, 0.0f, 1.0f));
      color[1]     = linearToSRGB(clamp(rgb[4*i+1]*a_normConst, 0.0f, 1.0f));
      color[2]     = linearToSRGB(clamp(rgb[4*i+2]*a_normConst, 0.0f, 1.0f));
      color[3]     = 1.0f;
      pixelData[i] = RealColorToUint32(color);
    }
  }
  else {
    const float invGamma  = 1.0f/a_gamma;
    for(int i=0;i<width*height;i++)
    {
      float color[4];
      color[0]     = clamp(std::pow(rgb[4*i+0]*a_normConst, invGamma), 0.0f, 1.0f);
      color[1]     = clamp(std::pow(rgb[4*i+1]*a_normConst, invGamma), 0.0f, 1.0f);
      color[2]     = clamp(std::pow(rgb[4*i+2]*a_normConst, invGamma), 0.0f, 1.0f);
      color[3]     = 1.0f;
      pixelData[i] = RealColorToUint32(color);
    }
  }
  return pixelData;
}


bool SaveImage4fToBMP(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, float a_gamma = 2.2f) 
{
  auto pixelData = FrameBufferColorToLDRImage(rgb,width,height,a_normConst,a_gamma);
  SaveBMP(outfilename, pixelData.data(), width, height);
  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> LoadImage4fFromEXR(const char* infilename, int* pW, int* pH) 
{
  std::vector<float> result;
  float* out; // width * height * RGBA
  int width  = 0;
  int height = 0;
  const char* err = nullptr; 

  int ret = LoadEXR(&out, &width, &height, infilename, &err);
  if (ret != TINYEXR_SUCCESS) {
    if (err) {
      fprintf(stderr, "[LoadImage4fFromEXR] : %s\n", err);
      std::cerr << "[LoadImage4fFromEXR] : " << err << std::endl;
      delete err;
    }
  }
  else {
    result.resize(width * height*4);
    *pW = uint32_t(width);
    *pH = uint32_t(height);
    memcpy(result.data(), out, width*height*sizeof(float)*4);
    free(out);
  }
  
  return result;
}

float* LoadImage4fFromEXRUnsafe(const char* infilename, int* pW, int* pH)
{
  float* out; // width * height * RGBA
  int width  = 0;
  int height = 0;
  const char* err = nullptr; 

  int ret = LoadEXR(&out, &width, &height, infilename, &err);
  if (ret != TINYEXR_SUCCESS) {
    if (err) {
      fprintf(stderr, "[LoadImage4fFromEXR] : %s\n", err);
      std::cerr << "[LoadImage4fFromEXR] : " << err << std::endl;
      delete err;
    }
    return nullptr;
  }
  else {
    *pW = uint32_t(width);
    *pH = uint32_t(height);
    return out;
  }
  
  return nullptr;
}

