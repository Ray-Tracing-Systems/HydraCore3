#ifndef INCLUDE_LITESCENE_SCENE_H_
#define INCLUDE_LITESCENE_SCENE_H_

#include <LiteScene/light_source.h>
#include <LiteScene/spectrum.h>
#include <LiteScene/texture.h>
#include <LiteScene/material.h>
#include <LiteScene/geometry.h>

#include <vector>
#include <string>

namespace ls {


    class HydraScene
    {
    public:
        std::vector<LightSource> light_sources;
        std::vector<Spectrum> spectra;
        std::vector<Texture> textures;
        std::vector<Material *> materials;
        std::vector<Geometry *> geometry;


        LightSource &get(SceneReference<LightSource> r) { return light_sources[r.id]; }
        const LightSource &get(SceneReference<LightSource> r) const { return light_sources[r.id]; }

        Spectrum &get(SceneReference<Spectrum> r) { return spectra[r.id]; }
        const Spectrum &get(SceneReference<Spectrum> r) const { return spectra[r.id]; }

        Texture &get(SceneReference<Texture> r) { return textures[r.id]; }
        const Texture &get(SceneReference<Texture> r) const { return textures[r.id]; }

        Material *get(SceneReference<Material> r) { return materials[r.id]; }
        const Material *get(SceneReference<Material> r) const { return materials[r.id]; }

        Geometry *get(SceneReference<Geometry> r) { return geometry[r.id]; }
        const Geometry *get(SceneReference<Geometry> r) const { return geometry[r.id]; }



        static uint32_t load(const std::string &path,       HydraScene &);
        static uint32_t save(const std::string &path, const HydraScene &);

    #ifdef LITESCENE_ENABLE_GLTF
        static uint32_t load_gltf(const std::string &path,       HydraScene &);
        static uint32_t save_gltf(const std::string &path, const HydraScene &);
    #endif

    };

}

#endif