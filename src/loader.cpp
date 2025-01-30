#include "loader.h"
#include "format.h"

#define VALIDATE_ID(x, name) \
    if((x) == SceneObject::INVALID_ID) { log_error("Id 0xFFFFFFFF is reserved (%s)", (name)); return ERROR_BAD_ID; }


namespace ls::internal {

    uint32_t SceneLoader::find_texture(const pugi::xml_node &node, std::optional<Texture *> &opt)
    {
        uint32_t id = node.attribute(L"id").as_uint();
        if(id != SceneObject::INVALID_ID) {
            auto it = textures.find(id);
            if(it != textures.end()) {
                opt = it->second;
                return SUCCESS;
            }
        }
        log_error("Unknown reference to texture %u", id);
        return ERROR_BAD_REFERENCE;
    }

    LiteMath::float4 SceneLoader::get_color(const pugi::xml_node& a_node)
    {
        auto val = hydra_xml::readvalVariant(a_node);
        if(std::holds_alternative<float>(val)) {
            return float4(std::get<float>(val));
        }
        else if(std::holds_alternative<float3>(val)) {
            return to_float4(std::get<float3>(val), 0.0f);
        }
        else {
            return std::get<float4>(val);
        }
    }

    uint32_t SceneLoader::find_spectrum(const pugi::xml_node& node, std::optional<Spectrum *> &opt)
    {
        uint32_t spec_id = node.attribute(L"id").as_uint();
        if(spec_id != SceneObject::INVALID_ID) {
            auto it = spectra.find(spec_id);
            if(it != spectra.end()) {
                opt = it->second;
                return SUCCESS;
            }
        }
        log_error("Invalid reference to spectrum %u", spec_id);
        return ERROR_BAD_REFERENCE;
    }

    uint32_t SceneLoader::load_color_holder(const pugi::xml_node &node, bool allow_spectrum, ColorHolder &clr)
    {
        clr.color = get_color(node);
        if(allow_spectrum) {
            auto specNode = node.child(L"spectrum");
            if(specNode) {
                return find_spectrum(specNode, clr.spectrum);
            }
        }
        return SUCCESS;
    }

    uint32_t SceneLoader::preload_lightsources(hydra_xml::HydraScene &scene)
    {
        uint32_t err;
        for(const auto &node : scene.LightNodes()) {
            const uint32_t id = node.attribute(L"id").as_uint();
            const std::string name = hydra_xml::ws2s(node.attribute(L"name").as_string());
            VALIDATE_ID(id, name.c_str());

            const std::wstring type = node.attribute(L"type").as_string();
            const std::wstring shape = node.attribute(L"shape").as_string();
            const std::wstring dist = node.attribute(L"distribution").as_string();

            LightSource *&lgt = light_sources[id];
            if(type == L"sky") {
                LightSourceSky *ptr = new LightSourceSky(name);
                const auto &texNode = node.child(L"intensity").child(L"color").child(L"texture");
                if(texNode) {
                    if(err = find_texture(texNode, ptr->texture)) {
                        return err;
                    }
                }
                const auto &backNode = node.child(L"back");
                if(backNode) {
                    if(err = find_texture(backNode, ptr->camera_back)) {
                        return err;
                    }
                }

                lgt = ptr;
            }
            else if(type == L"directional") {
                lgt = new LightSource(name, LightSourceType::DIRECTIONAL);
            }
            else if(shape == L"rect") {
                lgt = new LightSource(name, LightSourceType::RECT);
            }
            else if(shape == L"disk") {
                lgt = new LightSource(name, LightSourceType::DISK);
                lgt->radius = {node.child(L"size").attribute(L"radius").as_float()};
            }
            else if(shape == L"sphere") {
                lgt = new LightSource(name, LightSourceType::SPHERE);
                lgt->radius = {node.child(L"size").attribute(L"radius").as_float()};
            }
            else if(shape == L"point") {
                if(dist == L"spot") {
                    LightSourceSpot *ptr = new LightSourceSpot(name);
                    ptr->angle1 = hydra_xml::readval1f(node.child(L"falloff_angle"));
                    ptr->angle2 = hydra_xml::readval1f(node.child(L"falloff_angle2"));
                    const auto &projNode = node.child(L"projective");
                    if(projNode) {
                        SpotProj proj;
                        proj.fov = hydra_xml::readval1f(projNode.child(L"fov"));
                        proj.nearClipPlane = hydra_xml::readval1f(projNode.child(L"nearClipPlane"));
                        proj.farClipPlane = hydra_xml::readval1f(projNode.child(L"farClipPlane"));

                        const auto &projTexNode = projNode.child(L"texture");
                        if(projTexNode) {
                            if(err = find_texture(projTexNode, proj.texture)) {
                                return err;
                            }
                        }

                        ptr->projective = {proj};
                    }
                    lgt = ptr;
                }
                else {
                    lgt = new LightSource(name, LightSourceType::POINT);
                    if(dist == L"uniform" || dist == L"omni" || dist == L"ies") {
                        lgt->distribution = LightSourceDist::OMNI;
                    }
                }
            }
            else {
                log_error("Unknown light configuration : %s, %s (%s)", hydra_xml::ws2s(type).c_str(), hydra_xml::ws2s(shape).c_str(), name.c_str());
                return ERROR_LIGHTSOURCE_TYPE;
            }

            if(err = load_color_holder(node.child(L"intensity").child(L"color"), true, lgt->color)) {
                return err;
            }

            const auto &iesNode = node.child(L"ies");
            if(iesNode) {
                IES ies;
                ies.file_path = (scene_dir / hydra_xml::ws2s(iesNode.attribute(L"loc").as_string())).string();
                ies.point_area = iesNode.attribute(L"point_area").as_int() != 0;
                const auto &matrixAttrib = iesNode.attribute(L"matrix");
                if(matrixAttrib) {
                    ies.matrix = {hydra_xml::float4x4FromString(matrixAttrib.as_string())};
                }
                lgt->ies = {ies};
            }
        }

        return SUCCESS;
    }

    uint32_t SceneLoader::preload_geometry(hydra_xml::HydraScene &scene)
    {
        for(const auto &node : scene.GeomNodes()) {
            uint32_t id = node.attribute(L"id").as_uint();
            std::string name = hydra_xml::ws2s(node.attribute(L"name").as_string());

            VALIDATE_ID(id, name.c_str());

            std::wstring type = node.attribute(L"type").as_string();
            if(type == L"vsgf") {
                std::string loc = hydra_xml::ws2s(node.attribute(L"loc").as_string());
                geometry[id] = new Mesh(name, loc);
            }
            else {
                log_error("Unknown geometry type : %s", hydra_xml::ws2s(type).c_str());
                return ERROR_LIGHTSOURCE_TYPE;
            }
        }
        return SUCCESS;
    }


    uint32_t SceneLoader::preload(hydra_xml::HydraScene &scene)
    {
        uint32_t err;
        if(err = preload_geometry(scene)) return err;
        if(err = preload_spectra(scene)) return err;
        if(err = preload_textures(scene)) return err;
        if(err = preload_materials(scene)) return err;
        if(err = preload_lightsources(scene)) return err;
        return preload_instances(scene);
    }

}