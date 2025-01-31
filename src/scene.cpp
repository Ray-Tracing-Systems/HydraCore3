#include <LiteScene/scene.h>
#include "loader.h"
#include <filesystem>
#include <unordered_map>

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

        if(err = internal::move_map2vec(loader.geometry, scene.geometry)) return err;
        if(err = internal::move_map2vec(loader.spectra, scene.spectra)) return err;
        if(err = internal::move_map2vec(loader.textures, scene.textures)) return err;
        if(err = internal::move_map2vec(loader.materials, scene.materials)) return err;
        if(err = internal::move_map2vec(loader.light_sources, scene.light_sources)) return err;
        if(err = internal::move_map2vec(loader.scene_instances, scene.scene_instances, false)) return err;

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

