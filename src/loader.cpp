#include "loader.h"

namespace ls::internal {

    Texture *find_texture(Scene &data, const pugi::xml_node &node)
    {
        uint32_t id = node.attribute(L"id").as_uint();
        return data.textures[id];
    }

    ColorHolder load_color_holder(const pugi::xml_node &node, bool allow_spectrum)
    {

    }

    uint32_t preload_lightsources(Scene &data, const hydraxml::HydraScene &scene)
    {
        for(const auto &node : scene.LightNodes()) {
            const uint32_t id = node.attribute(L"id").as_uint();

            const std::wstring type = node.attribute(L"type").as_string();
            const std::wstring shape = node.attribute(L"shape").as_string();
            const std::wstring dist = node.attribute(L"distribution").as_string();

            LightSource *&lgt = data.light_sources[id];
            if(type == L"sky") {
                LightSourceSky *ptr = new LightSourceSky();
                const auto &texNode = node.child(L"intensity").child(L"color").child(L"texture");
                if(texNode) {
                    ptr->texture = {find_texture(data, texNode)};
                }
                const auto &backNode = node.child(L"back");
                if(backNode) {
                    ptr->camera_back = {find_texture(data, backNode)};
                }

                lgt = ptr;
            }
            else if(type == L"directional") {
                lgt = new LightSource(LightSourceType::DIRECTIONAL);
            }
            else if(shape == L"rect") {
                lgt = new LightSource(LightSourceType::RECT);
            }
            else if(shape == L"disk") {
                lgt = new LightSource(LightSourceType::DISK);
                lgt->radius = {node.child(L"size").attribute(L"radius").as_float()};
            }
            else if(shape == L"sphere") {
                lgt = new LightSource(LightSourceType::SPHERE);
                lgt->radius = {node.child(L"size").attribute(L"radius").as_float()};
            }
            else if(shape == L"point") {
                if(dist == L"spot") {
                    LightSourceSpot *ptr = new LightSourceSpot();
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
                            proj.texture = {find_texture(projTexNode)};
                        }

                        ptr->projective = {proj};
                    }
                    lgt = ptr;
                }
                else {
                    lgt = new LightSource(LightSourceType::POINT);
                    if(dist == L"uniform" || dist == L"omni" || dist == L"ies") {
                        lgt->distribution = LightSourceDist::OMNI;
                    }
                }
            }
            else {
                return ERROR_LIGHTSOURCE_TYPE;
            }

            lgt->color = load_color_holder(node.child(L"intensity").child(L"color"), true);


            const auto &iesNode = node.child(L"ies");
            if(iesNode) {
                IES ies;
                ies.file_path = (scene.scene_dir / hydra_xml::ws2s(iesNode.attribute(L"loc").as_string())).string();
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

}