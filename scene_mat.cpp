#include "scene.h"
#include "material.h"
#include "hydraxml.h"
#include "loadutil.h"
#include <memory>
#include <cstdio>


namespace LiteScene
{

    bool find_texture(const pugi::xml_node &colorNode, TextureInstance &inst)
    {
        auto texNode = colorNode.child(L"texture");
        if(texNode) {
            uint32_t id = texNode.attribute(L"id").as_uint();
            inst.id = id;
            return true;
        }
        return false;
    }

    void set_texture(pugi::xml_node &colorNode, const TextureInstance &inst)
    {

    }

    LiteMath::float4 get_color(const pugi::xml_node& a_node)
    {
        auto val = hydra_xml::readvalVariant(a_node);
        if(std::holds_alternative<float>(val)) {
            return LiteMath::float4(std::get<float>(val));
        }
        else if(std::holds_alternative<float3>(val)) {
            return LiteMath::to_float4(std::get<float3>(val), 0.0f);
        }
        else {
            return std::get<float4>(val);
        }
    }

    void set_color(pugi::xml_node &colorNode, const LiteMath::float4 &col)
    {
        if(col.w == 0.0f) {
            set_attr(colorNode, L"val", LM_to_wstring(LiteMath::float3(col.x, col.y, col.z)));
        }
        else {
            set_attr(colorNode, L"val", LM_to_wstring(col));
        }
    }

    bool load_color_holder(const pugi::xml_node &node, bool allow_spectrum, ColorHolder &clr)
    {
        clr.color = get_color(node);
        if(allow_spectrum) {
            auto specNode = node.child(L"spectrum");
            if(specNode) {
                uint32_t spec_id = specNode.attribute(L"id").as_uint();
                clr.spec_id = spec_id;
            }
        }
        return true;
    }

    void save_color_holder(pugi::xml_node &node, const ColorHolder &clr, bool allow_spectrum = true)
    {
        if(clr.color) {
            set_color(node, *clr.color);
        }
        if(clr.spec_id != INVALID_ID) {
            auto specNode = set_child(node, L"spectrum");
            set_attr(specNode, L"id", clr.spec_id);
            set_attr(specNode, L"type", L"ref");
        }
    }

    void load_gmc_node(const pugi::xml_node &node, std::variant<float, TextureInstance> &variant)
    {
        TextureInstance inst;
        if(find_texture(node, inst)) {
            variant = inst;
        }
        else {
            variant = node.attribute(L"val").as_float();
        }
    }

    Material *load_gltf_mat(uint32_t id, const std::string &name, const pugi::xml_node &node)
    {
        std::unique_ptr<GltfMaterial> mat{new GltfMaterial(id, name)};
        mat->raw_xml = node;

        //color
        auto colorNode = node.child(L"color");
        TextureInstance inst;
        if(find_texture(colorNode, inst)) {
            mat->color = std::move(inst);
        }
        else {
            auto col4 = get_color(colorNode);
            mat->color = LiteMath::float3(col4.x, col4.y, col4.z);
        }

        //gmc
        GltfMaterial::GMC gmc;
        bool validategmc = false;
        auto gmcNode = node.child(L"glossiness_metalness_coat");
        if(gmcNode) {
            TextureInstance gmcInst;
            if(find_texture(gmcNode, gmcInst)) {
                mat->glossiness_metalness_coat = gmcInst;
            }
            else {
                gmc.glossiness = gmc.metalness = gmc.coat = gmcNode.attribute(L"val").as_float();
                mat->glossiness_metalness_coat = gmc;
                validategmc = true;
            }
        }
        else {

            load_gmc_node(node.child(L"glossiness"), gmc.glossiness);
            load_gmc_node(node.child(L"metalness"), gmc.metalness);
            load_gmc_node(node.child(L"coat"), gmc.coat);
            mat->glossiness_metalness_coat = std::move(gmc);
            validategmc = true;
        }   

        //fresnel_ior
        mat->fresnel_ior = node.child(L"fresnel_ior").attribute(L"val").as_float();

        return mat.release();
    }

    Material *load_material(const pugi::xml_node &node)
    {
        uint32_t id = node.attribute(L"id").as_uint();
        std::string name = ws2s(node.attribute(L"name").as_string());

        const std::wstring type = node.attribute(L"type").as_string();
        if(type == L"gltf") {
            return load_gltf_mat(id, name, node);
        }
        else {
            printf("Material not supported\n");
            return nullptr;
        }

    }

    bool load_materials(HydraScene &scene, const pugi::xml_node &lib_node)
    {
        bool ok = true;

        for (pugi::xml_node node = lib_node.first_child(); node != nullptr && ok; node = node.next_sibling())
        {
            std::wstring type = node.attribute(L"type").as_string();


            if(std::wstring(node.name()) == L"material" && type == L"gltf")
            {   
                Material *mat;
                if(!(mat = load_material(node))) return false;
                scene.materials[mat->id] = mat;
            }
        }

        return ok;
    }

    void save_gmc_node(pugi::xml_node &node, const std::variant<float, TextureInstance> &val)
    {
        if(std::holds_alternative<TextureInstance>(val)) {
            set_texture(node, std::get<TextureInstance>(val));
        }
        else {
            set_attr(node, L"val", std::get<float>(val));
        }
    }

    bool save_gltf_mat(const GltfMaterial *mat, pugi::xml_node &node)
    {
        //color
        auto colorNode = set_child(node, L"color");
        if(std::holds_alternative<LiteMath::float3>(mat->color)) {
            set_color(colorNode, LiteMath::to_float4(std::get<LiteMath::float3>(mat->color), 0.0f));
        }
        else {
            set_texture(colorNode, std::get<TextureInstance>(mat->color));
        }

        //gmc
        if(std::holds_alternative<GltfMaterial::GMC>(mat->glossiness_metalness_coat)) {
            auto &gmc = std::get<GltfMaterial::GMC>(mat->glossiness_metalness_coat);
            auto glNode = set_child(node, L"glossiness");
            save_gmc_node(glNode, gmc.glossiness);
            auto meNode = set_child(node, L"metalness");
            save_gmc_node(meNode, gmc.metalness);
            auto coNode = set_child(node, L"coat");
            save_gmc_node(coNode, gmc.coat);
        }
        else {
            auto gmcNode = set_child(node, L"glossiness_metalness_coat");
            set_texture(gmcNode, std::get<TextureInstance>(mat->color));
        }

        auto fresNode = set_child(node, L"fresnel_ior");
        set_attr(fresNode, L"val", mat->fresnel_ior);
        return true;
    }

    bool save_material(const Material *mat, pugi::xml_node &node)
    {   
        set_attr(node, L"id", mat->id);
        set_attr(node, L"name", s2ws(mat->name));

        switch(mat->type()) {
        case MaterialType::GLTF:
            set_attr(node, L"type", L"gltf");
            return save_gltf_mat(static_cast<const GltfMaterial *>(mat), node);
        default:
            return false;
        }

        return true;
    }


    bool save_materials(const HydraScene &scene, pugi::xml_node lib_node)
    {
        for (const auto &[id, mat] : scene.materials)
        {
            auto node = lib_node.append_copy(mat->raw_xml);
            if(!save_material(mat, node)) return false;
        }
        return !lib_node.empty();
    }




}