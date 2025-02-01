#include "loader.h"
#include "format.h"

#define VALIDATE_ID(x, name) \
    if((x) == SceneObject::INVALID_ID) { log_error("Id 0xFFFFFFFF is reserved (%s)", (name)); return ERROR_BAD_ID; }


namespace ls::internal {


    static uint32_t load_old_hydra(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    } 

    static uint32_t load_lightsource_mat(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    }

    static uint32_t load_gltf_mat(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    }

    static uint32_t load_diffuse(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    }

    static uint32_t load_conductor(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    }

    static uint32_t load_dielectric(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    }

    static uint32_t load_plastic(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    }

    static uint32_t load_blend(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    }

    uint32_t load_thin_film(hydra_xml::HydraScene &scene, Material *&out)
    {
        return ERROR_NOT_IMPLEMENTED;
    }


    uint32_t SceneLoader::preload_materials(hydra_xml::HydraScene &scene)
    {
        uint32_t err;
        for(const auto &node : scene.MaterialNodes()) {
            uint32_t id = node.attribute(L"id").as_uint();
            const std::string name = hydra_xml::ws2s(node.attribute(L"name").as_string());
            VALIDATE_ID(id, name.c_str());

            const std::wstring type = node.attribute(L"type").as_string();
            if(type == L"hydra_material") {
                if(err = load_old_hydra(scene, materials[id])) return err;
            }
            else if(type == L"gltf") {
                if(err = load_gltf_mat(scene, materials[id])) return err;
            }
            else if(type == L"diffuse") {
                if(err = load_diffuse(scene, materials[id])) return err;
            }
            else if(type == L"rough_conductor") {
                if(err = load_conductor(scene, materials[id])) return err;
            }
            else if(type == L"dielectric") {
                if(err = load_dielectric(scene, materials[id])) return err;
            }
            else if(type == L"plastic") {
                if(err = load_plastic(scene, materials[id])) return err;
            }
            else if(type == L"blend") {
                if(err = load_blend(scene, materials[id])) return err;
            }
            else if(type == L"thin_film") {
                if(err = load_thin_film(scene, materials[id])) return err;
            }
            else {
                log_error("Unknown material type %s (%s)", hydra_xml::ws2s(type).c_str(), name.c_str());
                return ERROR_MATERIAL_TYPE;
            }
        }
        return SUCCESS;
    }

}