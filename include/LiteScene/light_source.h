#ifndef INCLUDE_LITESCENE_LIGHT__SOURCE_H_
#define INCLUDE_LITESCENE_LIGHT__SOURCE_H_
#include <LiteScene/sceneobj.h>
#include <LiteScene/material.h>
#include <optional>

namespace ls
{

    enum class LightSourceType
    {
        SKY, DIRECTIONAL,
        RECT, DISK, POINT, SPHERE //all in "area"
    };

    enum class LightSourceDist
    {
        LAMBERT, OMNI /* aka "omni", aka "ies" */, SPOT;
    };

    struct IES
    {
        std::string file_path;
        std::optional<LiteMath::float4x4> matrix;
        bool point_area;
    };

    class LightSource : public SceneObject
    {
    public:
       // bool visible;
        LightSourceDist distribution = LightSourceDist::LAMBERT;
        Material *material;
        ColorHolder color; //aka intensity.color
        float power; // aka intensity.multiplier
        std::optional<float> radius;
        std::optional<IES> ies;

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
        std::optional<Texture *> texture, camera_back;

        LightSourceSky()
            : LightSource(LightSourceType::SKY) { distribution = LightSourceDist::OMNI; }
    };

    struct SpotProj
    {
        float fov;
        float nearClipPlane;
        float farClipPlane;
        std::optional<Texture *> texture;
    };

    class LightSourceSpot : public LightSource
    {
    public:
        float angle1, angle2;
        std::optional<SpotProj> projective;

        LightSourceSpot()
            : LightSource(LightSourceType::POINT) { distribution = LightSourceDist::SPOT; }
    }

}

#endif