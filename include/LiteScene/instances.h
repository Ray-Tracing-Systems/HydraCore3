#ifndef INCLUDE_LITESCENE_INSTANCES_H_
#define INCLUDE_LITESCENE_INSTANCES_H_

#include <LiteScene/light_source.h>
#include <LiteScene/geometry.h>
#include <LiteMath.h>
#include <optional>
#include <vector>

namespace ls {

    class LightInstance : public SceneObject
    {
    public:
        LightSource *light;
        LiteMath::float4x4 matrix;
        uint32_t lgroup_id;
    };

    class SceneInstance;

    class GeometryInstance : public SceneObject
    {
    public:
        Geometry *object;
        uint32_t remap_id;
        uint32_t scene_sid;
        SceneInstance *scene_inst;
        LiteMath::float4x4 matrix;

        std::optional<LightInstance *> light_inst;
    };

    class SceneInstance : public SceneObject
    {
    public:
        std::vector<std::vector<uint32_t>> remaps;
        std::vector<LightInstance> light_instances;
        std::vector<GeometryInstance> geometry_instances;
    };


}

#endif