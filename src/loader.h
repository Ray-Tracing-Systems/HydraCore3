#ifndef SRC_LITESCENE_LOADER_H_
#define SRC_LITESCENE_LOADER_H_
#include <LiteScene/scene.h>
#include <map>
#include <filesystem>
#include "hydraxml.h"
#include <LiteMath.h>

namespace ls::internal {

    class SceneLoader
    {
    public:
        std::filesystem::path scene_dir;
        std::map<uint32_t, Geometry *> geometry;
        std::map<uint32_t, Spectrum *> spectra;
        std::map<uint32_t, Texture *> textures;
        std::map<uint32_t, Material *> materials;
        std::map<uint32_t, LightSource *> light_sources;
        std::map<uint32_t, SceneInstance *> scene_instances;

        uint32_t preload_lightsources(hydra_xml::HydraScene &scene);
        uint32_t preload_spectra(hydra_xml::HydraScene &scene) { return SUCCESS; }
        uint32_t preload_textures(hydra_xml::HydraScene &scene);
        uint32_t preload_materials(hydra_xml::HydraScene &scene);
        uint32_t preload_geometry(hydra_xml::HydraScene &scene);
        uint32_t preload_instances(hydra_xml::HydraScene &scene);

        uint32_t preload(hydra_xml::HydraScene &scene);

    private:
        LiteMath::float4 get_color(const pugi::xml_node& node);
        uint32_t find_texture(const pugi::xml_node &node, std::optional<Texture *> &); //textures must be preloaded
        uint32_t find_spectrum(const pugi::xml_node& node, std::optional<Spectrum *> &);
        uint32_t load_color_holder(const pugi::xml_node &node, bool allow_spectrum, ColorHolder &);
    };

    
    template<typename T>
    inline uint32_t move_map2vec(std::map<uint32_t, T *> &map, std::vector<T *> &vec, bool check_ids = true)
    {
        uint32_t vid = 0ul;
        for(const auto &[id, ptr] : map) {
            if(check_ids && id != vid) return ERROR_INCONSISTENT_ID;
            ptr->_set_id(vid++);
            vec.push_back(ptr);
        }
        map.clear();
        return SUCCESS;
    }

}


#endif