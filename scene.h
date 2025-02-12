#ifndef LITESCENE_SCENE_H_
#define LITESCENE_SCENE_H_
#include "LiteMath.h"
#include "3rd_party/pugixml.hpp"
#include "cmesh4.h"
#include "material.h"

#include <string>
#include <vector>
#include <map>

namespace LiteScene
{
    constexpr uint32_t    INVALID_ID = 0xFFFFFFFF;
    constexpr const char* INVALID_PATH = "INVALID_PATH";
    struct AABB
    {
        LiteMath::float3 boxMin;
        LiteMath::float3 boxMax;
    };

    // some context of the scene can affect how it's parts are stored/loaded
    // e.g. path to folder with .vsgf files is required to load geometry
    struct SceneMetadata
    {
        std::string scene_xml_folder = ""; //path to folder with .xml file
        std::string scene_xml_path = ""; //path to xml file. Empty if scene is created in memory and not saved yet
        std::string geometry_folder_relative = "./"; //path to folder with .vsgf files relative to scene_xml_path
        std::string geometry_folder = ""; //path to folder with .vsgf files, absolute or relative to current working directory

        pugi::xml_document xml_doc;
        pugi::xml_node custom_data; //all properties from scene xml that are not loaded to HydraScene
    };

    struct Spectrum 
    {

    };

    struct Texture
    {
        uint32_t id = INVALID_ID;
    };

    class Geometry
    {
    public:
        constexpr static uint32_t MESH_TYPE_ID   = 0;
        constexpr static uint32_t CUSTOM_TYPE_ID = 1000;

        virtual ~Geometry() {}
        bool load_node_base(pugi::xml_node node);
        pugi::xml_node save_node_base() const;
        virtual bool load_node(pugi::xml_node node) = 0;     // load properties from xml, no actual geometry data is loaded. Returns if loading was successful.
        virtual pugi::xml_node save_node() const = 0;        // save properties to xml.
        virtual bool load_data(const SceneMetadata &metadata) = 0;       // load actual data from file(s) (e.g. triangles from .vsgf file). Returns if loading was successful.
        virtual bool save_data(const SceneMetadata &metadata) = 0; // save actual data to file(s). Object should know how where to save. Returns if saving was successful.
        uint32_t id      = INVALID_ID;
        uint32_t type_id = INVALID_ID;
        uint32_t bytesize;
        std::string name;
        std::string type_name;

        pugi::xml_node custom_data; //all properties from xml node that are not loaded to struct fields
    };

    // classic mesh, not loaded by default, but can be loaded on request
    class MeshGeometry : public Geometry
    {
    public:
        virtual ~MeshGeometry() {}
        bool load_node(pugi::xml_node node) override;
        pugi::xml_node save_node() const override;
        bool load_data(const SceneMetadata &metadata) override;
        bool save_data(const SceneMetadata &metadata) override;

        bool is_loaded = false;
        std::string relative_file_path = INVALID_PATH;
        cmesh4::SimpleMesh mesh; // empty when not loaded
    };

    /* It is not a mesh, it is something else, like SDF or other implicit stuff
         By default, we just save its properties as a custom geometry and let user
         decide how to use it.*/
    class CustomGeometry : public Geometry
    {
    public:
        virtual ~CustomGeometry() {}
        bool load_node(pugi::xml_node node) override;
        pugi::xml_node save_node() const override;
        bool load_data(const SceneMetadata &metadata) override;
        bool save_data(const SceneMetadata &metadata) override;
    };

    struct LightSource
    {
        uint32_t id = INVALID_ID;
    };

    struct Camera
    {
        uint32_t id = INVALID_ID;
        std::string name;

        LiteMath::float3 pos;
        LiteMath::float3 lookAt;
        LiteMath::float3 up;
        float fov;
        float nearPlane;
        float farPlane;
        float exposureMult;
        LiteMath::float4x4 matrix;  // view matrix
        bool has_matrix;
        
        pugi::xml_node custom_data; //all properties from xml node that are not loaded to struct fields
    };

    struct RenderSettings
    {
        uint32_t id = INVALID_ID;
        std::string name;

        uint32_t width;
        uint32_t height;
        uint32_t spp;
        uint32_t depth;
        uint32_t depthDiffuse;
        pugi::xml_node custom_data; //all properties from xml node that are not loaded to struct fields
    };

    struct Instance
    {
        uint32_t id       = INVALID_ID;
        uint32_t mesh_id  = INVALID_ID;

        //idk if this is needed
        uint32_t rmap_id  = INVALID_ID;
        uint32_t scn_id   = INVALID_ID;
        uint32_t scn_sid  = INVALID_ID;
        uint32_t light_id = INVALID_ID;
        uint32_t linst_id = INVALID_ID;

        LiteMath::float4x4 matrix;                ///< transform matrix
        pugi::xml_node     custom_data; //all properties from xml node that are not loaded to struct fields
    };

    struct LightInstance
    {
        uint32_t id       = INVALID_ID;
        uint32_t mesh_id  = INVALID_ID;
        uint32_t light_id = INVALID_ID;

        //idk if this is needed
        uint32_t lgroup_id = INVALID_ID;

        LiteMath::float4x4 matrix;
        pugi::xml_node     custom_data; //all properties from xml node that are not loaded to struct fields
    };

    //InstancedScene is a weird concept, where there can be several scenes in one xml file
    //So each scene has it's own list of instances, bot lights and geometry.
    //Also remap lists are stored here
    struct InstancedScene
    {
        struct RemapList
        {
            uint32_t id = INVALID_ID;
            std::vector<uint32_t> remap;
        };

        uint32_t id = INVALID_ID;
        std::string name;
        AABB bbox;
        std::map<uint32_t, RemapList> remap_lists;
        std::map<uint32_t, Instance> instances;
        std::map<uint32_t, LightInstance> light_instances;
        pugi::xml_node     custom_data; //all properties from xml node that are not loaded to struct fields
    };

    struct HydraScene
    {
        HydraScene() {}
        ~HydraScene() { clear(); }
        //load scene from .xml file
        bool load(const std::string &filename); 

        //saves all the geometry to  a given folder and scene to xml file
        //it changes metadata, that's why it's not const
        bool save(const std::string &filename, const std::string &geometry_folder);

        //deletes all the data
        void clear();

        SceneMetadata metadata;

        std::map<uint32_t, Texture *> textures;  //HydraScene owns this data
        std::map<uint32_t, Material *> materials;//HydraScene owns this data
        std::map<uint32_t, Geometry *> geometries; //HydraScene owns this data
        std::map<uint32_t, LightSource> light_sources;
        std::map<uint32_t, Camera> cameras;
        std::map<uint32_t, RenderSettings> render_settings;
        std::map<uint32_t, InstancedScene> scenes;
    };
    std::wstring s2ws(const std::string& str);
    std::string ws2s(const std::wstring& wstr);
}

#endif