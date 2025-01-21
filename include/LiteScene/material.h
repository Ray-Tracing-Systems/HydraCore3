#ifndef INCLUDE_LITESCENE_MATERIAL_H_
#define INCLUDE_LITESCENE_MATERIAL_H_

#include <LiteScene/sceneobj.h>
#include <string>

namespace ls {


enum class MaterialType {
    GLTF,
    DIFFUSE,
    CONDUCTOR,
    PLASTIC,
    BLEND,
    THIN_FILM
};


class Material : public SceneObject
{
public:
    using SceneObject::SceneObject;


    virtual ~Material() = default;
    virtual MaterialType type() const = 0;
};


}

#endif