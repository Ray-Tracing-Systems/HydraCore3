#ifndef INCLUDE_LITESCENE_SCENE_H_
#define INCLUDE_LITESCENE_SCENE_H_
#include <LiteScene/material.h>
#include <LiteScene/mesh.h>
#include <vector>
#include <string>

namespace ls {


class HydraScene
{
public:

    static uint32_t load(const std::string &path,       HydraScene &);
    static uint32_t save(const std::string &path, const HydraScene &);
    static uint32_t load_gltf(const std::string &path,       HydraScene &);
    static uint32_t save_gltf(const std::string &path, const HydraScene &);

private:
    std::vector<SceneObject *> m_light_sources;
    std::vector<SceneObject *> m_spectra;
    std::vector<SceneObject *> m_textures;
    std::vector<Material *> m_materials;
    std::vector<Mesh> m_geometry;
};

}

#endif