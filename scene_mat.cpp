#include "scene.h"
#include "material.h"
#include "hydraxml.h"
#include "loadutil.h"
#include <memory>
#include <cstdio>


namespace LiteScene
{

    using LiteImage::Sampler;

    Sampler::AddressMode addr_mode_from_str(const std::wstring& a_mode)
    {
        if(a_mode == L"clamp")
            return Sampler::AddressMode::CLAMP;
        else if(a_mode == L"wrap")
            return Sampler::AddressMode::WRAP;
        else if(a_mode == L"mirror")
            return Sampler::AddressMode::MIRROR;
        else if(a_mode == L"border")
            return Sampler::AddressMode::BORDER;
        else if(a_mode == L"mirror_once")
            return Sampler::AddressMode::MIRROR_ONCE;
        else
            return Sampler::AddressMode::WRAP;
    }

    std::wstring addr_mode_to_str(Sampler::AddressMode mode)
    {
        switch(mode) {
        case Sampler::AddressMode::CLAMP:
            return L"clamp";
        case Sampler::AddressMode::MIRROR:
            return L"mirror";
        case Sampler::AddressMode::BORDER:
            return L"border";
        case Sampler::AddressMode::MIRROR_ONCE:
            return L"mirror_once";
        default:
            return L"wrap";
        }
    }

    bool load_texture_inst(const pugi::xml_node &texNode, TextureInstance &inst)
    {
        inst.id = texNode.attribute(L"id").as_uint();

        if(texNode.attribute(L"addressing_mode_u") != nullptr)
        {
            std::wstring addModeU = texNode.attribute(L"addressing_mode_u").as_string();
            inst.sampler.addr_mode_u  = addr_mode_from_str(addModeU);
        } 

        if(texNode.attribute(L"addressing_mode_v") != nullptr)
        {
            std::wstring addModeV = texNode.attribute(L"addressing_mode_v").as_string();
            inst.sampler.addr_mode_v  = addr_mode_from_str(addModeV);
        }

        if(texNode.attribute(L"addressing_mode_w") == nullptr)
            inst.sampler.addr_mode_w  = inst.sampler.addr_mode_v;
        else
        {
            std::wstring addModeW = texNode.attribute(L"addressing_mode_w").as_string();
            inst.sampler.addr_mode_w  = addr_mode_from_str(addModeW);
        }

        inst.sampler.filter = Sampler::Filter::LINEAR;
        if(texNode.attribute(L"filter") != nullptr)
        {
            std::wstring filterMode = texNode.attribute(L"filter").as_string();
            if(filterMode == L"point" || filterMode == L"nearest")
                inst.sampler.filter = Sampler::Filter::NEAREST;
            else if(filterMode == L"cubic" || filterMode == L"bicubic")
                inst.sampler.filter = Sampler::Filter::CUBIC;
        }

        if(texNode.attribute(L"input_gamma") != nullptr)
            inst.input_gamma = texNode.attribute(L"input_gamma").as_float();

        const std::wstring inputAlphaMode = texNode.attribute(L"input_alpha").as_string();
        if(inputAlphaMode == L"alpha") {
            inst.alpha_from_rgb = false;
        }

        inst.matrix = wstring_to_float4x4(texNode.attribute(L"matrix").as_string());
        return true;
    }

    bool find_texture(const pugi::xml_node &colorNode, TextureInstance &inst)
    {
        auto texNode = colorNode.child(L"texture");
        if(texNode) {
            return load_texture_inst(texNode, inst);
        }
        return false;
    }

    void set_texture(pugi::xml_node &colorNode, const TextureInstance &inst)
    {
        auto node = set_child(colorNode, L"texture");


        set_attr(node, L"id", inst.id);
        set_attr(node, L"type", L"texref");

        if(inst.sampler.addr_mode_u != Sampler::AddressMode::WRAP) {
            set_attr(node, L"addressing_mode_u", addr_mode_to_str(inst.sampler.addr_mode_u));
        }
        if(inst.sampler.addr_mode_v != Sampler::AddressMode::WRAP) {
            set_attr(node, L"addressing_mode_v", addr_mode_to_str(inst.sampler.addr_mode_v));
        }
        if(inst.sampler.addr_mode_w != Sampler::AddressMode::WRAP) {
            set_attr(node, L"addressing_mode_w", addr_mode_to_str(inst.sampler.addr_mode_w));
        }


        switch(inst.sampler.filter) {
        case Sampler::Filter::CUBIC:
            set_attr(node, L"filter", L"cubic");
            break;
        case Sampler::Filter::NEAREST:
            set_attr(node, L"filter", L"nearest");
            break;
        default:
            break;
        }

        if(inst.input_gamma != 2.2f) {
            set_attr(node, L"input_gamma", inst.input_gamma);
        }

        if(!inst.alpha_from_rgb) {
            set_attr(node, L"input_alpha", L"alpha");
        }

        set_attr(node, L"matrix", LM_to_wstring(inst.matrix));
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

    Material *convert_old_hydra(uint32_t id, const std::string &name, pugi::xml_node &node) {

        auto nodeEmiss = node.child(L"emission");
        if(nodeEmiss) { //Emissive
            std::unique_ptr<EmissiveMaterial> mat{new EmissiveMaterial(id, name)};

            

            return mat.release();
        }

        auto nodeDiffuse = node.child(L"diffuse");
        if(nodeDiffuse) { //Diffuse mat

        }

    }

    Material *load_material(pugi::xml_node &node)
    {
        uint32_t id = node.attribute(L"id").as_uint();
        std::string name = ws2s(node.attribute(L"name").as_string());

        const std::wstring type = node.attribute(L"type").as_string();
        if(type == L"gltf") {
            return load_gltf_mat(id, name, node);
        }
        else {
            Material *mat = new CustomMaterial(id, name);
            mat->raw_xml = node;
            return mat;
        }

    }

    bool load_materials(HydraScene &scene, pugi::xml_node &lib_node)
    {
        bool ok = true;

        for (pugi::xml_node node = lib_node.first_child(); node != nullptr && ok; node = node.next_sibling())
        {
            std::wstring type = node.attribute(L"type").as_string();


            if(std::wstring(node.name()) == L"material")
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
            return true;
        }

        return true;
    }


    bool save_materials(const HydraScene &scene, pugi::xml_node lib_node)
    {
        for (const auto &[id, mat] : scene.materials)
        {
            pugi::xml_node node = mat->raw_xml ? lib_node.append_copy(mat->raw_xml) : lib_node.append_child(L"material");
            if(!save_material(mat, node)) return false;
        }
        return !lib_node.empty();
    }




}