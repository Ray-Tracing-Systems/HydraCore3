#ifndef INCLUDE_LITESCENE_SCENEOBJ_H_
#define INCLUDE_LITESCENE_SCENEOBJ_H_
#include <cinttypes>
#include <string>
#include <vector>

namespace ls {

    class SceneObject
    {
    public:
        static inline uint32_t INVALID_ID = 0xFFFFFFFF;
        static inline uint32_t INLINE_ID = 0xFFFFFFFE;

        SceneObject(const std::string &name) : m_id(INVALID_ID), m_name(name) {}
        virtual ~SceneObject() = default;

        uint32_t id() const { return m_id; }
        std::string name() const { return m_name; }

        void _set_id(uint32_t id) { m_id = id; }
        
    private:
        uint32_t m_id;
        std::string m_name;
    };

}
#endif