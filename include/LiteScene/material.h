#ifndef INCLUDE_LITESCENE_MATERIAL_H_
#define INCLUDE_LITESCENE_MATERIAL_H_

#include <LiteScene/sceneobj.h>
#include <LiteScene/texture.h>
#include <LiteScene/spectrum.h>
#include <LiteMath.h>
#include <string>
#include <variant>
#include <optional>


namespace ls {

    enum class MaterialType 
    {
        GLTF,
        DIFFUSE,
        CONDUCTOR,
        DIELECTRIC,
        PLASTIC,
        BLEND,
        THIN_FILM
    };

    struct ColorHolder {
        std::optional<LiteMath::float4> color;
        std::optional<Spectrum *> spectrum;
    };

    class Material : public SceneObject
    {
    public:
        Material(const std::string &name) : SceneObject(name) {}

        virtual ~Material() = default;
        virtual MaterialType type() const = 0;
    };


    enum class AddressMode { CLAMP, MIRROR, BORDER, MIRROR_ONCE, WRAP };
    struct TextureInstance
    {
        SceneReference<Texture> texture;
        AddressMode addressing_mode_u;
        AddressMode addressing_mode_v;
        LiteMath::float4x4 matrix;
        //TODO ....
    };


    class GltfMaterial : public Material
    {
    public:

        std::variant<LiteMath::float3, TextureInstance> color;
        std::variant<float, TextureInstance> metalness;
        std::variant<float, TextureInstance> coat;
        std::variant<float, TextureInstance> glossiness;
        std::variant<float, TextureInstance> fresnel_ior;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::GLTF;
        }
    };

    class DiffuseMaterial : public Material
    {
    public:
        enum class BSDF { LAMBERT };

        std::variant<ColorHolder, TextureInstance> reflectance;
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
        std::variant<float, SceneReference<Spectrum>> eta;
        std::variant<float, SceneReference<Spectrum>> k;
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
        std::variant<float, SceneReference<Spectrum>> int_ior;
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
            std::variant<float, SceneReference<Spectrum>> eta;
            std::variant<float, SceneReference<Spectrum>> k;
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
        std::variant<float, SceneReference<Spectrum>> eta;
        std::variant<float, SceneReference<Spectrum>> k;
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
        SceneReference<Material> bsdf1;
        SceneReference<Material> bsdf2;

        using Material::Material;

        MaterialType type() const override
        {
            return MaterialType::BLEND;
        }
    };

}

#endif