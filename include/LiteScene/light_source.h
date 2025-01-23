#ifndef INCLUDE_LITESCENE_LIGHT__SOURCE_H_
#define INCLUDE_LITESCENE_LIGHT__SOURCE_H_
#include <LiteScene/sceneobj.h>
#include <LiteScene/material.h>

namespace ls
{

    enum class LightSourceType
    {
        SKY, DIRECTIONAL, RECT, DISK, POINT, SPHERE
    };

    enum class LightSourceDist
    {
        LAMBERT, OMNI, SPOT;
    };

    class LightSource : public SceneObject
    {
    public:
        using SceneObject::SceneObject;

        bool visible;
        LightSourceDist distribution;
        SceneReference<Material> material;


        virtual ~LightSource() = default;
        virtual LightSourceType type() const = 0;
    };

}

#endif