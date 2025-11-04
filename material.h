#ifndef LITESCENE_MATERIAL_H_
#define LITESCENE_MATERIAL_H_
#include "scene_common.h"
#include <string>
#include <variant>
#include <optional>


namespace LiteScene {

    struct TexSamplerHash
    {
        std::size_t operator()(const std::pair<TextureInstance::SamplerData, bool> &s) const 
        {
            std::size_t res = std::hash<uint64_t>{}(static_cast<uint64_t>(s.first.addr_mode_u)) * 31ull
                            + std::hash<uint64_t>{}(static_cast<uint64_t>(s.first.addr_mode_v)) * 31ull
                            + std::hash<uint64_t>{}(static_cast<uint64_t>(s.first.addr_mode_w)) * 31ull
                            + std::hash<uint64_t>{}(static_cast<uint64_t>(s.first.filter)) * 31ull
                            + std::hash<bool>{}(s.second);

            return res;
        }
    };



    enum class MaterialType 
    {
        CUSTOM,
        EMISSIVE,
        HYDRA_OLD,
        GLTF,
        DIFFUSE,
        CONDUCTOR,
        DIELECTRIC,
        PLASTIC,
        BLEND,
        THIN_FILM
    };

    class Material
    {
    public:
        Material(uint32_t _id, const std::string &_name) : id(_id), name(_name) {}

        virtual ~Material() = default;
        virtual MaterialType type() const = 0;

        uint32_t id;
        std::string name;
        pugi::xml_node raw_xml;
    };

    class CustomMaterial : public Material
    {
    public:
        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::CUSTOM;
        }
    };

    class EmissiveMaterial : public Material // Old Hydra for light sources
    {
    public:
        uint32_t light_id;
        std::variant<ColorHolder, TextureInstance> color;
        bool visible = true;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::EMISSIVE;
        }
    };

    class OldHydraMaterial : public Material // Old Hydra for other materials
    {
    public:
        using Material::Material;
        enum class BSDF { LAMBERT, OREN_NAYAR, NONE };
        enum class ReflBRDF { TORRANSE_SPARROW, PHONG, NONE };

        struct Data {
            LiteMath::float4 color;
            float glossiness;
            float ior;
        }; 

        std::variant<LiteMath::float4, TextureInstance> color = LiteMath::float4(0.0f, 0.0f, 0.0f, 0.0f);
        std::optional<float> diffuse_roughness;
        BSDF diffuse_bsdf_type = BSDF::NONE;
        ReflBRDF refl_brdf_type = ReflBRDF::NONE;

       
        std::optional<Data> reflectivity;
        std::optional<Data> transparency;


        MaterialType type() const override
        {
            return MaterialType::HYDRA_OLD;
        }
    };

    class GltfMaterial : public Material
    {
    public:

        struct GMC
        {
            std::variant<float, TextureInstance> glossiness;
            std::variant<float, TextureInstance> metalness;
            std::variant<float, TextureInstance> coat;
        };

        std::variant<LiteMath::float3, TextureInstance> color;
        std::variant<GMC, TextureInstance> glossiness_metalness_coat;
        float fresnel_ior;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::GLTF;
        }
    };

    class DiffuseMaterial : public Material
    {
    public:
        enum class BSDF { LAMBERT, OREN_NAYAR };

        std::variant<ColorHolder, TextureInstance> reflectance;
        float roughness = 0.0f; //Used only if BSDF::OPEN_NAYAR
        BSDF bsdf_type;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::DIFFUSE;
        }
    };

    struct MiRoughness
    {
        float alpha_u;
        float alpha_v;
    };

    class ConductorMaterial : public Material
    {
    public:
        enum class BSDF { GGX };
        BSDF bsdf_type;

        std::variant<MiRoughness, TextureInstance> alpha;
        std::variant<float, SceneRef<Spectrum>> eta;
        std::variant<float, SceneRef<Spectrum>> k;
        std::optional<LiteMath::float3> reflectance;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::CONDUCTOR;
        }
    };

    class DielectricMaterial : public Material
    {
    public:
        std::variant<float, SceneRef<Spectrum>> int_ior;
        float ext_ior;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::DIELECTRIC;
        }
    };

    class PlasticMaterial : public Material
    {
    public:
        std::variant<ColorHolder, TextureInstance> reflectance;
        float int_ior;
        float ext_ior;
        float alpha;
        bool nonlinear;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::PLASTIC;
        }
    };


    class ThinFilmMaterial : public Material
    {
    public:
        struct Layer
        {
            std::variant<float, SceneRef<Spectrum>> eta;
            std::variant<float, SceneRef<Spectrum>> k;
            float thickness;
        };
        struct ThicknessMap 
        {
            float min;
            float max;
            TextureInstance texture;
        };
        enum class BSDF { GGX };


        BSDF bsdf_type;
        std::variant<MiRoughness, TextureInstance> alpha;
        std::variant<float, SceneRef<Spectrum>> eta;
        std::variant<float, SceneRef<Spectrum>> k;
        float ext_ior;
        std::optional<ThicknessMap> thickness_map;
        bool transparent;
        std::vector<Layer> layers;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::THIN_FILM;
        }
    };

    class BlendMaterial : public Material
    {
    public:
        float weight;
        uint32_t bsdf1_id;
        uint32_t bsdf2_id;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::BLEND;
        }
    };

}

#endif