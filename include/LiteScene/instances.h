#ifndef INCLUDE_LITESCENE_INSTANCES_H_
#define INCLUDE_LITESCENE_INSTANCES_H_

#include <LiteScene/light_source.h>
#include <LiteScene/geometry.h>
#include <LiteMath.h>
#include <optional>
#include <vector>
#include <unordered_map>

namespace ls {

    class LightInstance : public SceneObject
    {
    public:
        LightSource *light;
        LiteMath::float4x4 matrix;
        uint32_t lgroup_id;

        LightInstance() : SceneObject("") {}
    };

    class SceneInstance;

    class GeometryInstance : public SceneObject
    {
    public:
        Geometry *object;
        std::optional<uint32_t> remap;
        //uint32_t scene_sid;
        //SceneInstance *scene_inst;
        LiteMath::float4x4 matrix;
        std::optional<LiteMath::float4x4> matrix_motion;

        std::optional<LightInstance *> light_inst;


        GeometryInstance() : SceneObject("") {}
    };

    class SceneInstance : public SceneObject
    {
    public:
        using SceneObject::SceneObject;

        std::vector<std::unordered_map<uint32_t, uint32_t>> remaps;
        std::vector<LightInstance *> light_instances;
        std::vector<GeometryInstance *> geometry_instances;

        ~SceneInstance()
        {
            for(auto *ptr : light_instances) delete ptr;
            for(auto *ptr : geometry_instances) delete ptr;
        }
    };


}

#endif