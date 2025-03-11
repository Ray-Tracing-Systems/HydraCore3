#include "scene.h"
#include "stb_image.h"
#include "stb_image_write.h"
#include <iostream>
#include <memory>

#define TINYGLTF_NO_INCLUDE_STB_IMAGE
#define TINYGLTF_NO_INCLUDE_STB_IMAGE_WRITE
#define TINYGLTF_IMPLEMENTATION
#include "3rd_party/tiny_gltf.h"

namespace gltf = tinygltf;



namespace LiteScene
{

    static void append_float3_as_float4(const gltf::Model &model, const gltf::Accessor& accessor, std::vector<LiteMath::float4> &target)
    {
        const gltf::BufferView& bufferView = model.bufferViews[accessor.bufferView];
        const gltf::Buffer& buffer = model.buffers[bufferView.buffer];
        const float* p = reinterpret_cast<const float*>(&buffer.data[bufferView.byteOffset + accessor.byteOffset]);
        for(size_t i = 0; i < accessor.count; ++i) {
            target.push_back(LiteMath::float4(
                                p[3 * i + 0],
                                p[3 * i + 1],
                                p[3 * i + 2],
                                0));
        }

    }

    template<typename T, typename P>
    void _append_int(const gltf::Model &model, const gltf::Accessor& accessor, std::vector<P> &target)
    {
        const gltf::BufferView& bufferView = model.bufferViews[accessor.bufferView];
        const gltf::Buffer& buffer = model.buffers[bufferView.buffer];

        const T* p = reinterpret_cast<const T*>(&buffer.data[bufferView.byteOffset + accessor.byteOffset]);
        target.insert(target.end(), p, p + accessor.count);
    }

    static void append_indices(const gltf::Model &model, const gltf::Accessor& accessor, std::vector<unsigned> &target)
    {
        switch(accessor.componentType) {
        case TINYGLTF_COMPONENT_TYPE_INT:
            _append_int<int32_t, unsigned>(model, accessor, target);
            return;
        case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
            _append_int<uint32_t, unsigned>(model, accessor, target);
            return;
        case TINYGLTF_COMPONENT_TYPE_SHORT:
            _append_int<int16_t, unsigned>(model, accessor, target);
            return;
        case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
            _append_int<uint16_t, unsigned>(model, accessor, target);
            return;
        case TINYGLTF_COMPONENT_TYPE_BYTE:
            _append_int<int8_t, unsigned>(model, accessor, target);
            return;
        case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
            _append_int<uint8_t, unsigned>(model, accessor, target);
            return;
        default:
            return;
        }

    }

    bool load_gltf_meshes(const std::string &filename, std::vector<Geometry *> out_meshes, std::vector<std::vector<uint32_t>> *materials)
    {
        gltf::Model model;
        gltf::TinyGLTF loader;
        std::string err, warn;

        if(!loader.LoadASCIIFromFile(&model, &err, &warn, filename)) {
            std::cerr << "Failed to parse glTF" << std::endl;
            return false;
        }

        if(!warn.empty()) {
            std::cout << "[Tiny-glTF WARN]: " << warn << std::endl;
        }
        if(!err.empty()) {
            std::cerr << "[Tiny-glTF ERROR]: " << err << std::endl;
            return false;
        }

        std::vector<std::unique_ptr<MeshGeometry>> meshes;
        meshes.reserve(model.meshes.size());


        for(const auto &mesh : model.meshes) {
            meshes.emplace_back(new MeshGeometry());
            MeshGeometry &mg = *(meshes.back());

            std::vector<uint32_t> prim_materials;

            cmesh4::SimpleMesh &simpleMesh = mg.mesh; 
            mg.is_loaded = true;

            for(const auto &prim : mesh.primitives) {
                if(prim.mode != TINYGLTF_MODE_TRIANGLES) {
                    std::cerr << "[ERROR] Only triangle primitives are supported" << std::endl;
                    return false;
                }

                const gltf::Accessor& posAccessor = model.accessors[prim.attributes.at("POSITION")];
                const gltf::Accessor& normAccessor = model.accessors[prim.attributes.at("NORMAL")];
                const gltf::Accessor& indAcessor  = model.accessors[prim.indices];
                append_float3_as_float4(model, posAccessor, simpleMesh.vPos4f);
                append_float3_as_float4(model, normAccessor, simpleMesh.vNorm4f);
                append_indices(model, indAcessor, simpleMesh.indices);

                if(materials) {
                    prim_materials.push_back(uint32_t(prim.material));
                }
            }
            if(materials) {
                materials->push_back(std::move(prim_materials));
            }
        }

        for(auto &mesh : meshes) out_meshes.push_back(mesh.release());

        return true;

    }

}