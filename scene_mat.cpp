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
            return LiteMath::to_float4(std::get<float3>(val), 1.0f);
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
            mat->raw_xml = node;

            mat->light_id = node.attribute(L"light_id").as_uint(); //???
            mat->visible = bool(node.attribute(L"visible").as_int());

            auto colorNode = nodeEmiss.child(L"color");
            TextureInstance inst;
            if(find_texture(colorNode, inst)) {
                mat->color = std::move(inst);
            }
            else {
                ColorHolder h;
                if(!load_color_holder(colorNode, true, h)) {
                    return nullptr;
                }
                mat->color = std::move(h);
            }

            return mat.release();
        }

        //Not light source

        std::unique_ptr<OldHydraMaterial> mat{new OldHydraMaterial(id, name)};
        mat->raw_xml = node;

        auto nodeDiffuse = node.child(L"diffuse");
        if(nodeDiffuse) { 
            std::wstring bsdf_type = nodeDiffuse.attribute(L"brdf_type"/*sic!*/).as_string();
            if(bsdf_type == L"orennayar") {
                mat->diffuse_bsdf_type = OldHydraMaterial::BSDF::OREN_NAYAR;
                auto nodeDiffRough = nodeDiffuse.child(L"roughness");
                if(nodeDiffRough) {
                    mat->diffuse_roughness = nodeDiffRough.attribute(L"val").as_float();
                }
            }
            else if(bsdf_type == L"lambert") {
                mat->diffuse_bsdf_type = OldHydraMaterial::BSDF::LAMBERT;
            }

            auto colorNode = nodeDiffuse.child(L"color");
            TextureInstance inst;
            if(find_texture(colorNode, inst)) {
                mat->color = std::move(inst);
            }
            else {
                mat->color = get_color(colorNode);
            }

        }

        auto nodeRefl = node.child(L"reflectivity");
        if(nodeRefl) {
            OldHydraMaterial::Data refl;
            refl.color = get_color(nodeRefl.child(L"color"));
            refl.glossiness = nodeRefl.child(L"glossiness").attribute(L"val").as_float();
            refl.ior = nodeRefl.child(L"fresnel_ior").attribute(L"val").as_float();
            mat->reflectivity = std::move(refl);

            std::wstring brdf_type = nodeRefl.attribute(L"brdf_type").as_string();
            if(brdf_type == L"torranse_sparrow") {
                mat->refl_brdf_type = OldHydraMaterial::ReflBRDF::TORRANSE_SPARROW;
            }
            else if(brdf_type == L"phong") {
                mat->refl_brdf_type = OldHydraMaterial::ReflBRDF::PHONG;
            }

        }

        auto nodeTransp = node.child(L"transparency");
        if(nodeTransp) {
            OldHydraMaterial::Data transp;
            transp.color = get_color(nodeTransp.child(L"color"));
            transp.glossiness = nodeTransp.child(L"glossiness").attribute(L"val").as_float();
            transp.ior = nodeTransp.child(L"ior").attribute(L"val").as_float();
            mat->transparency = std::move(transp);
        }
        return mat.release();
    }

    Material *load_diffuse_mat(uint32_t id, const std::string &name, pugi::xml_node &node)
    {
        std::unique_ptr<DiffuseMaterial> mat{new DiffuseMaterial(id, name)};
        const std::wstring bsdf_type = node.child(L"bsdf").attribute(L"type").as_string();

        if(bsdf_type == L"lambert") {
            mat->bsdf_type = DiffuseMaterial::BSDF::LAMBERT;
        }
        else if(bsdf_type == L"oren-nayar") {
            mat->bsdf_type = DiffuseMaterial::BSDF::OREN_NAYAR;
            auto nodeRoughness = node.child(L"roughness");
            if(nodeRoughness) {
                mat->roughness = nodeRoughness.attribute(L"val").as_float();
            }

        }
        else return nullptr;

        auto colorNode = node.child(L"reflectance");
        TextureInstance inst;
        if(find_texture(colorNode, inst)) {
            mat->reflectance = std::move(inst);
        }
        else {
            ColorHolder color;
            if(!load_color_holder(colorNode, true, color)) {
                return nullptr;
            }
            mat->reflectance = std::move(color);
        }

        return mat.release();

    }

    Material *load_material(pugi::xml_node &node)
    {
        uint32_t id = node.attribute(L"id").as_uint();
        std::string name = ws2s(node.attribute(L"name").as_string());

        const std::wstring type = node.attribute(L"type").as_string();
        if(type == L"gltf") {
            return load_gltf_mat(id, name, node);
        }
        if(type == L"hydra_material") {
            return convert_old_hydra(id, name, node);
        }
        if(type == L"diffuse") {
            return load_diffuse_mat(id, name, node);
        }
        
        Material *mat = new CustomMaterial(id, name);
        mat->raw_xml = node;
        return mat;


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

        set_val_child(node, L"fresnel_ior", mat->fresnel_ior);
        return true;
    }

    bool save_emissive_mat(const EmissiveMaterial *mat, pugi::xml_node &node)
    {
        set_attr(node, L"light_id", int(mat->light_id));
        set_attr(node, L"visible", int(mat->visible));

        auto nodeEmiss = set_child(node, L"emission");

        auto colorNode = set_child(nodeEmiss, L"color");
        if(std::holds_alternative<ColorHolder>(mat->color)) {
            save_color_holder(colorNode, std::get<ColorHolder>(mat->color), true);
        }
        else {
            set_texture(colorNode, std::get<TextureInstance>(mat->color));
        }

        return true;
    }

    bool has_default_color(const std::variant<LiteMath::float4, TextureInstance> &color)
    {
        return std::holds_alternative<LiteMath::float4>(color) && LiteMath::all_of(std::get<LiteMath::float4>(color) == LiteMath::float4(0, 0, 0, 0));
    }

    bool save_old_hydra_mat(const OldHydraMaterial *mat, pugi::xml_node &node)
    {   
        bool defColor = has_default_color(mat->color);
        if(mat->diffuse_bsdf_type != OldHydraMaterial::BSDF::NONE || !defColor) {
            auto nodeDiffuse = set_child(node, L"diffuse");
            if(mat->diffuse_bsdf_type == OldHydraMaterial::BSDF::LAMBERT) {
                set_attr(nodeDiffuse, L"brdf_type"/*sic!*/, L"lambert");   
            }
            else if(mat->diffuse_bsdf_type == OldHydraMaterial::BSDF::OREN_NAYAR) {
                set_attr(nodeDiffuse, L"brdf_type", L"orennayar");  
                if(mat->diffuse_roughness) {
                    set_val_child(nodeDiffuse, L"roughness", *mat->diffuse_roughness);
                }
            }

            if(!defColor) {
                auto colorNode = set_child(nodeDiffuse, L"color");
                if(std::holds_alternative<LiteMath::float4>(mat->color)) {
                    set_color(colorNode, std::get<LiteMath::float4>(mat->color));
                }
                else {
                    set_texture(colorNode, std::get<TextureInstance>(mat->color));
                }
            }
        }

        if(mat->reflectivity) {
            auto reflNode = set_child(node, L"reflectivity");
            const auto &refl = *mat->reflectivity;
            auto reflColorNode = set_child(reflNode, L"color");
            set_color(reflColorNode, refl.color);
            set_val_child(reflNode, L"glossiness", refl.glossiness);
            set_val_child(reflNode, L"fresnel_ior", refl.ior);

            if(mat->refl_brdf_type == OldHydraMaterial::ReflBRDF::TORRANSE_SPARROW) {
                set_attr(reflNode, L"brdf_type", L"torranse_sparrow");   
            }
            else if(mat->refl_brdf_type == OldHydraMaterial::ReflBRDF::PHONG) {
                set_attr(reflNode, L"brdf_type", L"phong");  
            }

        }

        if(mat->transparency) {
            auto transpNode = set_child(node, L"transparency");
            const auto &transp = *mat->transparency;
            auto transpColorNode = set_child(transpNode, L"color");
            set_color(transpColorNode, transp.color);
            set_val_child(transpNode, L"glossiness", transp.glossiness);
            set_val_child(transpNode, L"ior", transp.ior);
        }

        return true;
    }

    bool save_diffuse_mat(const DiffuseMaterial *mat, pugi::xml_node &node)
    {
        auto bsdfNode = set_child(node, L"bsdf");
        if(mat->bsdf_type == DiffuseMaterial::BSDF::LAMBERT) {
            set_attr(bsdfNode, L"type", "lambert");
        }
        else {
            set_attr(bsdfNode, L"type", "oren-nayar");
            if(mat->roughness != 0.0f) {
                set_val_child(node, L"roughness", mat->roughness);
            }
        }


        auto colorNode = set_child(node, L"reflectance");
        if(std::holds_alternative<TextureInstance>(mat->reflectance)) {
            set_texture(colorNode, std::get<TextureInstance>(mat->reflectance));
        }
        else {
            save_color_holder(colorNode, std::get<ColorHolder>(mat->reflectance));
        }

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
        case MaterialType::EMISSIVE:
            set_attr(node, L"type", L"hydra_material");
            return save_emissive_mat(static_cast<const EmissiveMaterial *>(mat), node);
        case MaterialType::HYDRA_OLD:
            set_attr(node, L"type", L"hydra_material");
            return save_old_hydra_mat(static_cast<const OldHydraMaterial *>(mat), node);
        case MaterialType::DIFFUSE:
            set_attr(node, L"type", L"diffuse");
            return save_diffuse_mat(static_cast<const DiffuseMaterial *>(mat), node);
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