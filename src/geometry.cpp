#include <LiteScene/geometry.h>
#include <fstream>

namespace {

    cmesh4::Header header(const std::string &path)
    {
        std::ifstream str(path);
        return cmesh4::LoadHeader(str);
    }

}

namespace ls {

    Mesh::Mesh(const std::string &name, std::shared_ptr<cmesh4::SimpleMesh> mesh) : Geometry(name), m_mesh(std::move(mesh)), m_header() {}

    Mesh::Mesh(const std::string &name, const std::string &path) : Geometry(name), m_mesh(path), m_header(header(path)) {}

    std::shared_ptr<cmesh4::SimpleMesh const> Mesh::get_mesh() const
    {
        if(std::holds_alternative<std::shared_ptr<cmesh4::SimpleMesh>>(m_mesh)) {
            auto &ptr = std::get<std::shared_ptr<cmesh4::SimpleMesh>>(m_mesh);
            return std::const_pointer_cast<cmesh4::SimpleMesh const>(ptr);
        }
        else {
            const std::string &path = std::get<std::string>(m_mesh);
            cmesh4::SimpleMesh mesh = cmesh4::LoadMeshFromVSGF(path.c_str());
            return std::shared_ptr<cmesh4::SimpleMesh const>(new cmesh4::SimpleMesh(std::move(mesh)));
        }
    }

    std::shared_ptr<cmesh4::SimpleMesh> Mesh::get_mesh()
    {
        if(std::holds_alternative<std::shared_ptr<cmesh4::SimpleMesh>>(m_mesh)) {
            return std::get<std::shared_ptr<cmesh4::SimpleMesh>>(m_mesh);
        }
        else {
            const std::string &path = std::get<std::string>(m_mesh);
            cmesh4::SimpleMesh mesh = cmesh4::LoadMeshFromVSGF(path.c_str());
            return std::shared_ptr<cmesh4::SimpleMesh>(new cmesh4::SimpleMesh(std::move(mesh)));
        }
    }

    void Mesh::set_mesh(std::shared_ptr<cmesh4::SimpleMesh> new_mesh)
    {
        m_header = {};
        m_mesh = std::move(new_mesh);
    }


    void Mesh::set_mesh(const std::string &file)
    {
        m_header = header(file);
        m_mesh = file;
    }


}