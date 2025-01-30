#ifndef INCLUDE_LITESCENE_SCENE_H_
#define INCLUDE_LITESCENE_SCENE_H_

#include <LiteScene/light_source.h>
#include <LiteScene/spectrum.h>
#include <LiteScene/texture.h>
#include <LiteScene/material.h>
#include <LiteScene/geometry.h>
#include <LiteScene/instances.h>
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

        const std::vector<LightSource *> &get_lightsources() { return light_sources; }
        const std::vector<Spectrum *> &get_spectra() { return spectra; }
        const std::vector<Texture *> &get_textures() { return textures; }
        const std::vector<Material *> &get_materials() { return materials; }
        const std::vector<Geometry *> &get_geometrys() { return geometry; }

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
    constexpr uint32_t ERROR_NOT_IMPLEMENTED = 0xFFFFFFFF;
    constexpr uint32_t ERROR_XML_LOADER = 1;
    constexpr uint32_t ERROR_BAD_ID = 2;
    constexpr uint32_t ERROR_BAD_REFERENCE = 3;
    constexpr uint32_t ERROR_LIGHTSOURCE_TYPE = 4;
    constexpr uint32_t ERROR_GEOMETRY_TYPE = 5;
    constexpr uint32_t ERROR_MATERIAL_TYPE = 6;






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