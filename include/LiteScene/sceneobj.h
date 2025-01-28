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

        SceneObject() : m_id(INVALID_ID), m_name() {}
        virtual ~SceneObject() = default;

        uint32_t id() const { return m_id; }
        std::wstring name() const { return m_name; }
        
    private:
        uint32_t m_id;
        std::wstring m_name;

        friend class HydraScene;
    };

}
#endif