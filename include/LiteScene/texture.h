#ifndef INCLUDE_LITESCENE_TEXTURE_H_
#define INCLUDE_LITESCENE_TEXTURE_H_

#include <LiteScene/sceneobj.h>


namespace ls {

class Texture : public SceneObject
{
public:
    using SceneObject::SceneObject;

    std::shared_ptr<LiteImage::ICombinedImageSampler> sampler;
};

enum class AddressMode { CLAMP, WRAP, MIRROR, BORDER, MIRROR_ONCE, WRAP; }

struct TextureInstance
{
    SceneReference<Texture> texture;
    AddressMode addressing_mode_u;
    AddressMode addressing_mode_v;
    
};

#endif