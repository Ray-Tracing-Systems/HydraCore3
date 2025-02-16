#include "scene.h"
#include "loadutil.h"

#include <filesystem>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#ifndef LITESCENE_NOTINYEXR
#define TINYEXR_IMPLEMENTATION
#endif
#include "tinyexr.h"

namespace fs = std::filesystem;
using LiteImage::Image2D;
using namespace LiteMath;

namespace LiteScene
{
    std::vector<float> LoadImage1fFromEXR(const char* infilename, int* pW, int* pH)
    {
        float* out; // width * height * RGBA
        int width       = 0;
        int height      = 0;
        const char* err = nullptr;

        int ret = LoadEXR(&out, &width, &height, infilename, &err);
        if (ret != TINYEXR_SUCCESS) {
            if (err) {
                fprintf(stderr, "[LoadImage1fFromEXR] : %s\n", err);
                std::cerr << "[LoadImage1fFromEXR] : " << err;
                std::cerr << " from path : " << infilename << std::endl;
                
                delete err;
            }
            return std::vector<float>();
        }

        const int imgSize = width * height;
        std::vector<float> result(imgSize);
        *pW = uint32_t(width);
        *pH = uint32_t(height);
        
        #pragma omp parallel for
        for(int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                size_t idx = (x + (height - y - 1) * width) * 4;
                size_t out_idx = x + y * width;
                if (std::isinf(out[idx]))
                    result[out_idx] = 65504.0f;                       // max half float according to ieee
                else
                    result[out_idx] = clamp(out[idx], 0.0f, 65504.0f); // max half float according to ieee
                
            }
        }

        free(out);
        return result;
    }


    std::vector<float> LoadImage4fFromEXR(const char* infilename, int* pW, int* pH) 
    {
        float* out; // width * height * RGBA
        int width       = 0;
        int height      = 0;
        const char* err = nullptr; 

        int ret = LoadEXR(&out, &width, &height, infilename, &err);
        if (ret != TINYEXR_SUCCESS) {
            if (err) {
                fprintf(stderr, "[LoadImage4fFromEXR] : %s\n", err);
                std::cerr << "[LoadImage4fFromEXR] : " << err << std::endl;
                delete err;
            }
            return std::vector<float>();
        }

        std::vector<float> result(width * height * 4);
        *pW = uint32_t(width);
        *pH = uint32_t(height);

        #pragma omp parallel for
        for(int y = 0; y < height; y++)
        {
            const int offset1 = (height - y - 1) * width * 4;
            const int offset2 = y * width * 4;
            memcpy((void*)(result.data() + offset1), (void*)(out + offset2), width * sizeof(float) * 4);
        }
        free(out);  
        
        return result;
    }


    float* LoadImage4fFromEXRUnsafe(const char* infilename, int* pW, int* pH)
    {
        float* out; // width * height * RGBA
        int width  = 0;
        int height = 0;
        const char* err = nullptr; 

        int ret = LoadEXR(&out, &width, &height, infilename, &err);
        if (ret != TINYEXR_SUCCESS) {
            if (err) {
                fprintf(stderr, "[LoadImage4fFromEXR] : %s\n", err);
                std::cerr << "[LoadImage4fFromEXR] : " << err << std::endl;
                delete err;
            }
            return nullptr;
        }

        *pW = uint32_t(width);
        *pH = uint32_t(height);
        return out;
    }

    std::shared_ptr<LiteImage::ICombinedImageSampler> make_combined_sampler(Texture::Info &info, const TextureInstance &inst)
    {
        std::shared_ptr<LiteImage::ICombinedImageSampler> pResult;

        const bool disable_gamma = inst.input_gamma == 1.0f;
        LiteImage::Sampler sampler {
            .addressU = inst.sampler.addr_mode_u,
            .addressV = inst.sampler.addr_mode_v,
            .addressW = inst.sampler.addr_mode_w,
            .filter = inst.sampler.filter
        };
        const fs::path path{info.path};
        const std::string ex = path.extension().string();

        if(ex == ".bmp" || ex == ".ppm" || ex == ".jpg" || ex == ".jpeg" || ex == ".png") {

            Image2D<uint32_t> image = LiteImage::LoadImage<uint32_t>(info.path.c_str());
            auto pTexture = std::make_shared<Image2D<uint32_t>>(std::move(image));
            pTexture->setSRGB(!disable_gamma);
            pResult = LiteImage::MakeCombinedTexture2D(pTexture, sampler);

        }
        else if (ex == ".exr") {
            int wh[2]{ 0, 0 };
            if (info.bpp == 16) {
                const auto image_vect = LoadImage4fFromEXR(info.path.c_str(), &wh[0], &wh[1]);
                auto pTexture = std::make_shared<Image2D<float4>>(wh[0], wh[1], (const float4*)image_vect.data());
                pTexture->setSRGB(false);
                pResult = LiteImage::MakeCombinedTexture2D(pTexture, sampler);
            }
            else {
                const auto image_vect = LoadImage1fFromEXR(info.path.c_str(), &wh[0], &wh[1]);
                auto pTexture = std::make_shared<Image2D<float>>(wh[0], wh[1], (const float*)image_vect.data());
                pTexture->setSRGB(false);
                pResult = LiteImage::MakeCombinedTexture2D(pTexture, sampler);
            }
        }
        else if(ex.find(".image") != std::string::npos) { // hydra image formats: image4f, image4ub

            int wh[2] = {0,0};

            std::ifstream fin(info.path.c_str(), std::ios::binary);
            if(!fin.is_open())
              std::cout << "[LoadTextureAndMakeCombined]: can't open '" << info.path << "'" << std::endl;
            
            fin.read((char*)wh, sizeof(int)*2);
            if(wh[0] == 0 || wh[1] == 0) {
                float4 data[1] = {float4(1.0f, 1.0f, 1.0f, 1.0f)};
                auto pTexture  = std::make_shared<Image2D<float4>>(1, 1, data);
                pTexture->setSRGB(false);
                pResult = LiteImage::MakeCombinedTexture2D(pTexture, sampler);
            }
            else if(info.bpp == 16) { // image4f
                std::vector<float> data(4 * wh[0] * wh[1]);
                fin.read((char*)data.data(), 4 * sizeof(float) * data.size());
                fin.close();

                auto pTexture = std::make_shared<Image2D<float4>>(wh[0], wh[1], (const float4*)data.data());
                pResult = LiteImage::MakeCombinedTexture2D(pTexture, sampler);
            }
            else {                       // image4ub
                std::vector<uint32_t> data(wh[0]*wh[1]);
                fin.read((char*)data.data(), sizeof(uint32_t)*data.size());
                fin.close();

                auto pTexture = std::make_shared< Image2D<uint32_t> >(wh[0], wh[1], data.data());
                pTexture->setSRGB(!disable_gamma);
                pResult = LiteImage::MakeCombinedTexture2D(pTexture, sampler);
            }
        }

        return pResult;
    }

    std::shared_ptr<LiteImage::ICombinedImageSampler> Texture::get_combined_sampler(const TextureInstance &inst)
    {
        std::pair<TextureInstance::SamplerData, bool> sd = {inst.sampler, inst.input_gamma == 1.0f};
        auto it = tex_cache.find(sd);
        if(it != tex_cache.end()) {
            return it->second;
        }
        auto s = make_combined_sampler(info, inst);
        tex_cache.insert({sd, s});
        return s;
    }

    bool Texture::load_info(pugi::xml_node &node, const std::string &scene_root)
    {
        if(node.attribute(L"loc").empty()) {
            info.path = ws2s(std::wstring(node.attribute(L"path").as_string()));
        }
        else {
            std::string loc = ws2s(std::wstring(node.attribute(L"loc").as_string()));
            info.path = (fs::path(scene_root) / loc).string();
        }
        
        info.width = node.attribute(L"width").as_uint();
        info.height = node.attribute(L"height").as_uint();
        if(info.width != 0 && info.height != 0) {
            const size_t byteSize = node.attribute(L"bytesize").as_ullong();
            info.bpp = uint32_t(byteSize / size_t(info.width * info.height));
        }
        return true;
    }

    //https://stackoverflow.com/questions/60530422/how-to-check-a-file-is-contained-in-a-folder-with-c
    // =======
    fs::path normalized_trimed(const fs::path& p)
    {
        auto r = p.lexically_normal();
        if (r.has_filename()) return r;
        return r.parent_path();
    }

    fs::path get_relative_if_possible(const fs::path& dir, const fs::path& target)
    {   
        if(!dir.empty()) {
            auto base = normalized_trimed(dir);
            auto sub = normalized_trimed(target).parent_path();
            auto m = std::mismatch(base.begin(), base.end(), 
                                   sub.begin(), sub.end());

            if(m.first == base.end()) { // dir contains target
                fs::path rel_path = std::accumulate(m.second, sub.end(), fs::path{}, std::divides{}) / target.filename();
                return rel_path.lexically_normal();
            }
        }
        return target.lexically_normal();
    }
    // =======

    bool Texture::save_info(pugi::xml_node &node, const std::string &scene_root) const
    {
        fs::path path = get_relative_if_possible(fs::path(scene_root), fs::path(info.path));
        if(path.is_absolute()) {
            set_attr(node, L"path", s2ws(path));
        }
        else {
            set_attr(node, L"loc", s2ws(path));
        }

        set_attr(node, L"height", info.height);
        set_attr(node, L"width", info.width);
        if(info.width != 0 && info.height != 0) {
            set_attr(node, L"bytesize", size_t(info.bpp) * info.width * info.height);
        }
        return true;
    }

}