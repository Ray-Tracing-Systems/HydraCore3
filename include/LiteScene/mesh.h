#ifndef INCLUDE_LITESCENE_MESH_H_
#define INCLUDE_LITESCENE_MESH_H_
#include <LiteScene/sceneobj.h>
#include <LiteScene/cmesh4.h>

namespace ls {

class Mesh : public SceneObject
{
public:
    const cmesh4::SimpleMesh &get_mesh() const { return m_mesh; }
    const LiteMath::float4x4 &get_transform(uint32_t id) const { return m_transforms[id]; }
private:
    cmesh4::SimpleMesh m_mesh;
    std::vector<LiteMath::float4x4> m_transforms;
};

}

#endif