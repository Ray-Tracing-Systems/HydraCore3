#include <LiteScene/scene.h>
#include "loader.h"
#include <filesystem>

namespace {

    template<typename T>
    void move_map2vec(std::unordered_map<uint32_t, T *> &map, std::vector<T *> &vec)
    {
        uint32_t vid = 0ul;
        for(const auto &[id, ptr] : map) {
            ptr->_set_id(vid++);
            vec.push_back(ptr);
        }
        map.clear();
    }

}

namespace ls {

    uint32_t HydraScene::load(const std::string &scenePath, HydraScene &scene)
    {
        uint32_t err;
        hydra_xml::HydraScene scene_xml;
        if(scene_xml.LoadState(scenePath)) return ERROR_XML_LOADER;
        
        std::filesystem::path path{scenePath};
        internal::SceneLoader loader;
        loader.scene_dir = path.parent_path().string();
        if(err = loader.preload(scene_xml)) {
            return err;
        }

        move_map2vec(loader.geometry, scene.geometry);
        move_map2vec(loader.spectra, scene.spectra);
        move_map2vec(loader.textures, scene.textures);
        move_map2vec(loader.materials, scene.materials);
        move_map2vec(loader.light_sources, scene.light_sources);
        move_map2vec(loader.scene_instances, scene.scene_instances);
        
        return SUCCESS;
    }



    uint32_t HydraScene::save(const std::string &path, const HydraScene &)
    {
        return 0;
    }

    HydraScene::~HydraScene()
    {
        for(auto *ptr : light_sources) delete ptr;
        for(auto *ptr : spectra) delete ptr;
        for(auto *ptr : textures) delete ptr;
        for(auto *ptr : materials) delete ptr;
        for(auto *ptr : geometry) delete ptr;
        for(auto *ptr : scene_instances) delete ptr;
    }

}

