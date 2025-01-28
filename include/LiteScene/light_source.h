#ifndef INCLUDE_LITESCENE_LIGHT__SOURCE_H_
#define INCLUDE_LITESCENE_LIGHT__SOURCE_H_
#include <LiteScene/sceneobj.h>
#include <LiteScene/material.h>

namespace ls
{

    enum class LightSourceType
    {
        SKY, DIRECTIONAL,
        RECT, DISK, POINT, SPHERE //all in "area"
    };

    enum class LightSourceDist
    {
        LAMBERT, UNIFORM /* aka "omni", aka "ies" */, SPOT;
    };

    class LightSource : public SceneObject
    {
    public:
        bool visible;
        LightSourceDist distribution;
        SceneReference<Material> material;
        ColorHolder color; //aka intensity.color
        float power; // aka intensity.multiplier

        LightSource(LightSourceType type)
            : m_type(type) {}


        virtual ~LightSource() = default;
        LightSourceType type() const { return m_type; }
    private:
        LightSourceType m_type;
    };

    class LightSourceSky : public LightSource
    {
    public:
        std::optional<SceneReference<Texture>> texture, camera_back;

        LightSourceSky()
            : LightSource(LightSourceType::SKY) {}
    };

}

#endif