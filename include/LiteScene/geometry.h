#ifndef INCLUDE_LITESCENE_GEOMETRY_H_
#define INCLUDE_LITESCENE_GEOMETRY_H_

#include <LiteScene/sceneobj.h>
#include <string>
#include <LiteScene/cmesh4.h>


namespace ls {


enum class GeometryType {
    MESH
};


class Geometry : public SceneObject
{
public:
    using SceneObject::SceneObject;


    virtual ~Geometry() = default;
    virtual GeometryType type() const = 0;
};


class Mesh : public Geometry
{
public:
    const cmesh4::SimpleMesh &get_mesh() const { return m_mesh; }
    const LiteMath::float4x4 &get_transform(uint32_t id) const { return m_transforms[id]; }

    GeometryType type() const override { return GeometryType::MESH; }

private:
    cmesh4::SimpleMesh m_mesh;
    std::vector<LiteMath::float4x4> m_transforms;
};

#endif