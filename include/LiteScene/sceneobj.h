#ifndef INCLUDE_LITESCENE_SCENEOBJ_H_
#define INCLUDE_LITESCENE_SCENEOBJ_H_
#include <cinttypes>
#include <string>
#include <vector>

namespace ls {

    template<typename T>
    struct SceneReference
    {
        uint32_t id;
    };


    class SceneObject
    {
    public:
        static inline uint32_t INVALID_ID = 0xFFFFFFFF;
        static inline uint32_t INLINE_ID = 0xFFFFFFFE;

        SceneObject(const std::string &name) : m_id(INVALID_ID), m_name(name) {}
        virtual ~SceneObject() = default;

        uint32_t id() const { return m_id; }
        std::string name() const { return m_name; }
        
    private:
        uint32_t m_id;
        std::string m_name;

        friend class HydraScene;
    };

}
#endif