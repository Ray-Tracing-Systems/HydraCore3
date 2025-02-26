#ifndef LITESCENE_SCENE__COMMON_H_
#define LITESCENE_SCENE__COMMON_H_
#include "3rd_party/pugixml.hpp"
#include "LiteMath.h"
#include "Image2d.h"
#include <optional>


namespace LiteScene
{
    constexpr uint32_t INVALID_ID = 0xFFFFFFFF;
    constexpr char* INVALID_PATH = "INVALID_PATH";

    struct ColorHolder {
        std::optional<LiteMath::float4> color;
        uint32_t spec_id = INVALID_ID;
    };

    struct TextureInstance
    {
        struct SamplerData {
            LiteImage::Sampler::AddressMode addr_mode_u;
            LiteImage::Sampler::AddressMode addr_mode_v;
            LiteImage::Sampler::AddressMode addr_mode_w;
            LiteImage::Sampler::Filter filter;

            bool operator==(const SamplerData &other) const
            {
                return addr_mode_u == other.addr_mode_u
                    && addr_mode_v == other.addr_mode_v
                    && addr_mode_w == other.addr_mode_w
                    && filter == other.filter;
            }
        };

        uint32_t id;

        SamplerData sampler;
        LiteMath::float4x4 matrix;
        float input_gamma = 2.2f;
        bool alpha_from_rgb = true;
    };

    template<typename T>
    struct SceneRef {
        uint32_t id;
        operator uint32_t() const { return id; }
    };




    class Texture;
    class Material;
    class Spectrum;
    class Geometry;
    class LightSource;
}

#endif