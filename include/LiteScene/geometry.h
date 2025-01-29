#ifndef INCLUDE_LITESCENE_GEOMETRY_H_
#define INCLUDE_LITESCENE_GEOMETRY_H_

#include <LiteScene/sceneobj.h>
#include <LiteScene/cmesh4.h>
#include <variant>
#include <memory>
#include <optional>


namespace ls {

    enum class GeometryType 
    {
        MESH
    };


    class Geometry : public SceneObject
    {
    public:

        virtual ~Geometry() = default;
        virtual GeometryType type() const = 0;
    };


    class Mesh : public Geometry
    {
    public:
        Mesh(std::shared_ptr<cmesh4::SimpleMesh> mesh);
        Mesh(const std::string &path);

        std::shared_ptr<cmesh4::SimpleMesh const> get_mesh() const;
        void set_mesh(std::shared_ptr<cmesh4::SimpleMesh> new_mesh);
        void set_mesh(const std::string &file);

        const std::optional<cmesh4::Header> &metadata_if_file() const { return m_header; }

        GeometryType type() const override { return GeometryType::MESH; }

    private:
        std::variant<
            std::shared_ptr<cmesh4::SimpleMesh>,
            std::string
        > m_mesh;

        std::optional<cmesh4::Header> m_header;
    };

}

#endif