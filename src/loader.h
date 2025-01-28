#ifndef SRC_LITESCENE_LOADER_H_
#define SRC_LITESCENE_LOADER_H_
#include <LiteScene/scene.h>
#include <unordered_map>
#include <filesystem>
#include "hydraxml.h"

namespace ls::internal {

    struct Scene
    {
        std::filesystem::path scene_dir;
        std::unordered_map<uint32_t, LightSource *> light_sources;
        std::unordered_map<uint32_t, Spectrum *> spectra;
        std::unordered_map<uint32_t, Texture *> textures;
        std::unordered_map<uint32_t, Material *> materials;
        std::unordered_map<uint32_t, Geometry *> geometry;
        std::unordered_map<uint32_t, SceneInstance *> scene_instances;
    };

    ColorHolder load_color_holder(const pugi::xml_node &node, bool allow_spectrum);
    Texture *find_texture(const pugi::xml_node &node); //textures must be preloaded

    uint32_t preload_lightsources(Scene &data, const hydraxml::HydraScene &scene);
    uint32_t preload_spectra(Scene &data, const hydraxml::HydraScene &scene);
    uint32_t preload_textures(Scene &data, const hydraxml::HydraScene &scene);
    uint32_t preload_materials(Scene &data, const hydraxml::HydraScene &scene);
    uint32_t preload_geometry(Scene &data, const hydraxml::HydraScene &scene);
    uint32_t preload_instances(Scene &data, const hydraxml::HydraScene &scene);
}


#endif