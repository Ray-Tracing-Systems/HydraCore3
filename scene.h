#ifndef INCLUDE_LITESCENE_SCENE_H_
#define INCLUDE_LITESCENE_SCENE_H_

#include <LiteScene/light_source.h>
#include <LiteScene/spectrum.h>
#include <LiteScene/texture.h>
#include <LiteScene/material.h>
#include <LiteScene/geometry.h>
#include <LiteScene/instances.h>
#include <unordered_map>
#include <vector>
#include <memory>
#include <string>

namespace ls {

    class HydraScene
    {
    public:

        template<typename T, typename ...Args>
        T *make_lightsource(Args &&...args);
        template<typename T, typename ...Args>
        T *make_spectrum(Args &&...args);
        template<typename T, typename ...Args>
        T *make_texture(Args &&...args);
        template<typename T, typename ...Args>
        T *make_material(Args &&...args);
        template<typename T, typename ...Args>
        T *make_geometry(Args &&...args);

        LightSource *get_lightsource(uint32_t id) { return light_sources[id]; }
        const LightSource *get_lightsource(uint32_t id) const { return light_sources[id]; }

        Spectrum *get_spectrum(uint32_t id) { return spectra[id]; }
        const Spectrum *get_spectrum(uint32_t id) const { return spectra[id]; }

        Texture *get_texture(uint32_t id) { return textures[id]; }
        const Texture *get_texture(uint32_t id) const { return textures[id]; }

        Material *get_material(uint32_t id) { return materials[id]; }
        const Material *get_material(uint32_t id) const { return materials[id]; }

        Geometry *get_geometry(uint32_t id) { return geometry[id]; }
        const Geometry *get_geometry(uint32_t id) const { return geometry[id]; }

        //NOTE: currently only one scene instance is supported
        SceneInstance *get_scene(uint32_t id = 0) { return scene_instances[id]; }
        const SceneInstance *get_scene(uint32_t id = 0) const { return scene_instances[id]; }

        const std::vector<LightSource *> &get_lightsources() { return light_sources; }
        const std::vector<Spectrum *> &get_spectra() { return spectra; }
        const std::vector<Texture *> &get_textures() { return textures; }
        const std::vector<Material *> &get_materials() { return materials; }
        const std::vector<Geometry *> &get_geometries() { return geometry; }
        const std::vector<SceneInstance *> &get_scene_instances() { return scene_instances; }

        static uint32_t load(const std::string &path,       HydraScene &);
        static uint32_t save(const std::string &path, const HydraScene &);

#ifdef LITESCENE_ENABLE_GLTF
        static uint32_t load_gltf(const std::string &path,       HydraScene &);
        static uint32_t save_gltf(const std::string &path, const HydraScene &);
#endif
        ~HydraScene();

    private:
        std::vector<LightSource *> light_sources;
        std::vector<Spectrum *> spectra;
        std::vector<Texture *> textures;
        std::vector<Material *> materials;
        std::vector<Geometry *> geometry;
        std::vector<SceneInstance *> scene_instances;
    };

    constexpr uint32_t SUCCESS = 0;
    constexpr uint32_t ERROR_NOT_IMPLEMENTED = 0;//0xFFFFFFFF;
    constexpr uint32_t ERROR_XML_LOADER = 1;
    constexpr uint32_t ERROR_BAD_ID = 2;
    constexpr uint32_t ERROR_BAD_REFERENCE = 3;
    constexpr uint32_t ERROR_INCONSISTENT_ID = 4;

    constexpr uint32_t ERROR_GEOMETRY_TYPE = 100;
    constexpr uint32_t ERROR_MATERIAL_TYPE = 200;
    constexpr uint32_t ERROR_LIGHTSOURCE_TYPE = 300;
    
    constexpr uint32_t ERROR_BAD_RMAPLIST = 601;



    template<typename T, typename ...Args>
    T *HydraScene::make_lightsource(Args &&...args)
    {
        T *obj = new T(std::forward<Args...>(args...));
        light_sources.push_back(obj);
        obj->_set_id(light_sources.size() - 1);
        return obj;
    }

    template<typename T, typename ...Args>
    T *HydraScene::make_spectrum(Args &&...args)
    {
        T *obj = new T(std::forward<Args...>(args...));
        spectra.push_back(obj);
        obj->_set_id(spectra.size() - 1);
        return obj;
    }

    template<typename T, typename ...Args>
    T *HydraScene::make_texture(Args &&...args)
    {
        T *obj = new T(std::forward<Args...>(args...));
        textures.push_back(obj);
        obj->_set_id(textures.size() - 1);
        return obj;
    }

    template<typename T, typename ...Args>
    T *HydraScene::make_material(Args &&...args)
    {
        T *obj = new T(std::forward<Args...>(args...));
        materials.push_back(obj);
        obj->_set_id(materials.size() - 1);
        return obj;
    }

    template<typename T, typename ...Args>
    T *HydraScene::make_geometry(Args &&...args)
    {
        T *obj = new T(std::forward<Args...>(args...));
        geometry.push_back(obj);
        obj->_set_id(geometry.size() - 1);
        return obj;
    }

}

#endif