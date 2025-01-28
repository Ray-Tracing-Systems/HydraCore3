#include <LiteScene/scene.h>
#include "hydraxml.h"

namespace {
    uint32_t load_light_source()
}

namespace ls {

    uint32_t HydraScene::load(const std::string &path, HydraScene &scene)
    {
        hydra_xml::HydraScene scene_xml;
        if(scene_xml.LoadState(path)) return ERROR_XML_LOADER;
        
        

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

