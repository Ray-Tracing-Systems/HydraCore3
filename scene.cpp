#include "scene.h"
#include "hydraxml.h"
#include "loadutil.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <cassert>
#include <locale>
#include <codecvt>

namespace LiteScene
{

    std::wstring s2ws(const std::string& str)
    {
        using convert_typeX = std::codecvt_utf8<wchar_t>;
        std::wstring_convert<convert_typeX, wchar_t> converterX;
        return converterX.from_bytes(str);
    }

    std::string ws2s(const std::wstring& wstr)
    {
        using convert_typeX = std::codecvt_utf8<wchar_t>;
        std::wstring_convert<convert_typeX, wchar_t> converterX;
        return converterX.to_bytes(wstr);
    }

    std::vector<std::string> split(std::string aStr, char aDelim)
    {
        std::vector<std::string> res;
        std::string str = aStr.substr(0, aStr.find(aDelim));

        while (str.size() < aStr.size())
        {
            res.push_back(str);
            aStr = aStr.substr(aStr.find(aDelim) + 1);
            str = aStr.substr(0, aStr.find(aDelim));
        }

        res.push_back(str);

        return res;
    }

    std::wstring save_float_array_to_string(const std::vector<float> &array)
    {
        std::wstringstream outputStream;
        for (int i = 0; i < array.size(); i++)
        {
            outputStream << array[i] << L" ";
        }
        return outputStream.str();
    }

  
    std::wstring float4x4ToString(const LiteMath::float4x4 &matrix)
    {
        return LM_to_wstring(matrix);
        /*return save_float_array_to_string({ matrix.get_row(0).x, matrix.get_row(0).y, matrix.get_row(0).z, matrix.get_row(0).w,
                                                                             matrix.get_row(1).x, matrix.get_row(1).y, matrix.get_row(1).z, matrix.get_row(1).w,
                                                                             matrix.get_row(2).x, matrix.get_row(2).y, matrix.get_row(2).z, matrix.get_row(2).w,
                                                                             matrix.get_row(3).x, matrix.get_row(3).y, matrix.get_row(3).z, matrix.get_row(3).w });
    */}

    AABB AABBFromString(const std::wstring &str)
    {
        auto data = wstring_to_float_arr(str, 6);
        AABB result;
        result.boxMin = LiteMath::float3(data[0], data[1], data[2]);
        result.boxMax = LiteMath::float3(data[3], data[4], data[5]);
        return result;
    }

    std::wstring AABBToString(const AABB &aabb)
    {
        return save_float_array_to_string({ aabb.boxMin.x, aabb.boxMin.y, aabb.boxMin.z, aabb.boxMax.x, aabb.boxMax.y, aabb.boxMax.z });
    }

    bool Geometry::load_node_base(pugi::xml_node node)
    {
        custom_data = node;

        id = node.attribute(L"id").as_uint(INVALID_ID);
        bytesize = node.attribute(L"bytesize").as_uint(0);
        name = ws2s(node.attribute(L"name").as_string());
        type_name = ws2s(node.attribute(L"type").as_string());

        if (id == INVALID_ID)
        {
            printf("[Geometry::load_node_base] Invalid id\n");
            return false;
        }

        return true;
    }
    
    void Geometry::save_node_base(pugi::xml_node &node) const
    {        
        set_attr(node, L"id", id);
        set_attr(node, L"bytesize", bytesize);
        set_attr(node, L"name", s2ws(name));
        set_attr(node, L"type", s2ws(type_name));
    }



    bool MeshGeometry::load_node(pugi::xml_node node)
    {
        bool ok = load_node_base(node);
        if (!ok) return false;

        if (node.attribute(L"loc").empty())
        {
            relative_file_path = INVALID_PATH;
            printf("[MeshGeometry::load_node] No location\n");
            return false;
        }
        relative_file_path = ws2s(node.attribute(L"loc").as_string());
        type_id = MESH_TYPE_ID;

        return true;
    }

    bool MeshGeometry::save_node(pugi::xml_node &node) const
    {
        if (relative_file_path == INVALID_PATH)
        {
            printf("[MeshGeometry::save_node] No location is specified. Save data first\n");
            return false;
        }
        save_node_base(node);
        node.set_name(L"mesh");
        set_attr(node, L"loc", s2ws(relative_file_path));
        return true;
    }
    bool MeshGeometry::load_data(const SceneMetadata &metadata)
    {
        if (is_loaded)
            return true;

        if (relative_file_path == INVALID_PATH)
        {
            printf("[MeshGeometry::load_data] No location is specified. Load node first\n");
            return false;
        }
        std::string path = metadata.scene_xml_folder + "/" + relative_file_path;
        mesh = cmesh4::LoadMeshFromVSGF(path.c_str());
        bool ok = mesh.VerticesNum() > 0 && mesh.IndicesNum() > 0;
        if (!ok)
        {
            printf("[MeshGeometry::load_data] Failed to load mesh %s\n", path.c_str());
            return false;
        }
        is_loaded = true;
        return true;
    }

    bool MeshGeometry::save_data(const SceneMetadata &metadata)
    {
        if (!is_loaded)
        {
            printf("[MeshGeometry::save_data] Mesh is not loaded\n");
            return false;
        }
        std::string name = "mesh_" + std::to_string(id);
        relative_file_path = metadata.geometry_folder_relative + "/" + name + ".vsgf";
        std::string file_path = metadata.scene_xml_folder == "" ? relative_file_path : 
        metadata.scene_xml_folder + "/" + relative_file_path;
        cmesh4::SaveMeshToVSGF(file_path.c_str(), mesh);

        return true;
    }

    bool CustomGeometry::load_node(pugi::xml_node node)
    {
        bool ok = load_node_base(node);
        if (!ok) return false;

        type_id = CUSTOM_TYPE_ID;

        return true;    
    }

    bool CustomGeometry::save_node(pugi::xml_node &node) const
    {
        save_node_base(node);
        return true;
    }
    bool CustomGeometry::load_data(const SceneMetadata &metadata)
    {
        printf("[CustomGeometry::load_data] Cannot load data for custom geometry of unknown type\n");
        return false;
    }

    bool CustomGeometry::save_data(const SceneMetadata &metadata)
    {
        printf("[CustomGeometry::save_data] Cannot save data for custom geometry of unknown type\n");
        return false;
    }

    void HydraScene::initialize_empty_scene()
    {
        metadata.xml_doc.load_string(LR""""(
            <?xml version="1.0"?>
            <textures_lib />
            <materials_lib />
            <geometry_lib />
            <lights_lib />
            <cam_lib />
            <render_lib />
            <scenes />
        )"""");
    }

    void HydraScene::clear()
    {
        //for (int i = 0; i < textures.size(); i++)
        //    delete textures[i];

        for (int i = 0; i < materials.size(); i++)
            delete materials[i];

        for (int i = 0; i < geometries.size(); i++)
            delete geometries[i];
        
        for (int i = 0; i < light_sources.size(); i++)
            delete light_sources[i];


        metadata = SceneMetadata();

        textures.clear();
        materials.clear();
        geometries.clear();

        light_sources.clear();
        cameras.clear();
        render_settings.clear();
        scenes.clear();

        initialize_empty_scene();
    }

    bool load_geometry(HydraScene &scene, pugi::xml_node lib_node)
    {
        bool ok = true;

        for (pugi::xml_node geom_node = lib_node.first_child(); geom_node != nullptr && ok; geom_node = geom_node.next_sibling())
        {
            if (std::wstring(geom_node.name()) == L"mesh")
            {
                MeshGeometry *geom = new MeshGeometry();
                ok = geom->load_node(geom_node);
                if (ok)
                    scene.geometries[geom->id] = geom;
            }
            else
            {
                CustomGeometry *geom = new CustomGeometry();
                ok = geom->load_node(geom_node);
                if (ok)
                    scene.geometries[geom->id] = geom;
            }
        }

        if (scene.geometries.size() == 0)
        {
            printf("[HydraScene::load_geometry] No geometries loaded\n");
            return false;
        }

        return true;
    }

    bool load_camera(Camera &cam, const pugi::xml_node &cam_node)
    {
        cam.id = cam_node.attribute(L"id").as_uint();
        cam.name = ws2s(cam_node.attribute(L"name").as_string());

        cam.fov       = hydra_xml::readval1f(cam_node.child(L"fov")); 
        cam.nearPlane = hydra_xml::readval1f(cam_node.child(L"nearClipPlane"));
        cam.farPlane  = hydra_xml::readval1f(cam_node.child(L"farClipPlane"));  

        auto expNode = cam_node.child(L"exposure_mult");
        if(expNode)
            cam.exposureMult = hydra_xml::readval1f(expNode);  
        else
            cam.exposureMult = 1.0f;
        
        cam.pos    = hydra_xml::readval3f(cam_node.child(L"position"));
        cam.lookAt = hydra_xml::readval3f(cam_node.child(L"look_at"));
        cam.up     = hydra_xml::readval3f(cam_node.child(L"up"));

        cam.has_matrix = false;
        if(cam_node.child(L"matrix"))
        {
            cam.matrix = LiteMath::transpose(wstring_to_float4x4(cam_node.child(L"matrix").attribute(L"val").as_string()));
            cam.has_matrix = true;
        }

        cam.custom_data = cam_node;

        return true;
    }

    bool load_cameras(HydraScene &scene, const pugi::xml_node &lib_node)
    {
        bool ok = true;

        for (pugi::xml_node cam_node = lib_node.first_child(); cam_node != nullptr && ok; cam_node = cam_node.next_sibling())
        {
            if(std::wstring(cam_node.name()) == L"camera")
            {
                Camera camera;
                if(!load_camera(camera, cam_node)) return false;
                scene.cameras[camera.id] = std::move(camera);
            }

        }

        return ok;
    }

    bool load_render_settings(RenderSettings &settings, const pugi::xml_node &node)
    {
        settings.id = node.attribute(L"id").as_uint();
        settings.name = ws2s(node.attribute(L"name").as_string());

        settings.width  = hydra_xml::readval1u(node.child(L"width"));
        settings.height = hydra_xml::readval1u(node.child(L"height"));
        settings.depth  = hydra_xml::readval1u(node.child(L"trace_depth"));
        settings.depthDiffuse = hydra_xml::readval1u(node.child(L"diff_trace_depth"));
        settings.spp    = hydra_xml::readval1u(node.child(L"maxRaysPerPixel"));
        settings.custom_data = node;
        return true;
    }

    bool load_all_render_settings(HydraScene &scene, const pugi::xml_node &lib_node)
    {
        bool ok = true;

        for (pugi::xml_node set_node = lib_node.first_child(); set_node != nullptr && ok; set_node = set_node.next_sibling())
        {
            if(std::wstring(set_node.name()) == L"render_settings")
            {
                RenderSettings settings;
                if(!load_render_settings(settings, set_node)) return false;
                scene.render_settings[settings.id] = std::move(settings);
            }

        }

        return ok;
    }

    bool load_textures(HydraScene &scene, const pugi::xml_node &lib_node)
    {
        bool ok = true;

        for (pugi::xml_node tex_node = lib_node.first_child(); tex_node != nullptr && ok; tex_node = tex_node.next_sibling())
        {
            if(std::wstring(tex_node.name()) == L"texture")
            {
                Texture tex;
                if(!tex.load_info(tex_node, scene.metadata.scene_xml_folder)) return false;
                scene.textures[tex.id] = std::move(tex);
            }

        }

        return ok;
    }

    bool find_texture(const pugi::xml_node &colorNode, TextureInstance &inst);
    bool load_color_holder(const pugi::xml_node &node, bool allow_spectrum, ColorHolder &clr);

    LightSource *load_lightsource(pugi::xml_node &node, const SceneMetadata &meta)
    {
        const uint32_t id = node.attribute(L"id").as_uint();
        const std::string name = ws2s(node.attribute(L"name").as_string());
        const uint32_t mat_id = node.attribute(L"mat_id").as_uint();

        const std::wstring type = node.attribute(L"type").as_string();
        const std::wstring shape = node.attribute(L"shape").as_string();
        const std::wstring dist = node.attribute(L"distribution").as_string();

        std::unique_ptr<LightSource> lgt;
        if(type == L"sky") {
            LightSourceSky *ptr = new LightSourceSky();
            const auto &colNode = node.child(L"intensity").child(L"color");
            if(colNode) {
                TextureInstance inst;
                if(find_texture(colNode, inst)) {
                    ptr->texture = std::move(inst);
                }
            }
            const auto &backNode = node.child(L"back");
            if(backNode) {
                TextureInstance inst;
                if(find_texture(backNode, inst)) {
                    ptr->camera_back = std::move(inst);
                }
            }

            lgt.reset(ptr);
        }
        else if(type == L"directional") {
            lgt.reset(new LightSource(LightSource::Type::DIRECTIONAL));
        }
        else if(shape == L"rect") {
            lgt.reset(new LightSource(LightSource::Type::RECT));
            auto sizeNode = node.child(L"size");
            lgt->half_width = sizeNode.attribute(L"half_width").as_float();
            lgt->half_length = sizeNode.attribute(L"half_length").as_float();
        }
        else if(shape == L"disk") {
            lgt.reset(new LightSource(LightSource::Type::DISK));
            lgt->radius = node.child(L"size").attribute(L"radius").as_float();
        }
        else if(shape == L"sphere") {
            lgt.reset(new LightSource(LightSource::Type::SPHERE));
            lgt->radius = node.child(L"size").attribute(L"radius").as_float();
        }
        else if(shape == L"point") {
            if(dist == L"spot") {
                LightSourceSpot *ptr = new LightSourceSpot();
                ptr->angle1 = hydra_xml::readval1f(node.child(L"falloff_angle"));
                ptr->angle2 = hydra_xml::readval1f(node.child(L"falloff_angle2"));
                const auto &projNode = node.child(L"projective");
                if(projNode) {
                    LightSourceSpot::Proj proj;
                    proj.fov = hydra_xml::readval1f(projNode.child(L"fov"));
                    proj.nearClipPlane = hydra_xml::readval1f(projNode.child(L"nearClipPlane"));
                    proj.farClipPlane = hydra_xml::readval1f(projNode.child(L"farClipPlane"));

                    const auto &projTexNode = projNode.child(L"texture");
                    if(projTexNode) {
                        TextureInstance inst;
                        if(find_texture(projNode, inst)) {
                            proj.texture = std::move(inst);
                        }
                    }

                    ptr->projective = std::move(proj);
                }
                lgt.reset(ptr);
            }
            else {
                lgt.reset(new LightSource(LightSource::Type::POINT));
                if(dist == L"uniform" || dist == L"omni" || dist == L"ies") {
                    lgt->distribution = LightSource::Dist::OMNI;
                }
            }
        }
        else {
            return nullptr;
        }

        auto intNode = node.child(L"intensity");
        if(!load_color_holder(intNode.child(L"color"), true, lgt->color)) {
            return nullptr;
        }
        lgt->power = intNode.child(L"multiplier").attribute(L"val").as_float();

        const auto &iesNode = node.child(L"ies");
        if(iesNode) {
            LightSource::IES ies;
            ies.file_path = (fs::path(meta.scene_xml_folder) / ws2s(iesNode.attribute(L"loc").as_string())).string();
            ies.point_area = iesNode.attribute(L"point_area").as_int() != 0;
            const auto &matrixAttrib = iesNode.attribute(L"matrix");
            if(matrixAttrib) {
                ies.matrix = wstring_to_float4x4(matrixAttrib.as_string());
            }
            lgt->ies = std::move(ies);
        }

        lgt->id = id;
        lgt->name = std::move(name);
        lgt->mat_id = mat_id;
        return lgt.release();
    }

    bool load_lightsources(HydraScene &scene, const pugi::xml_node &lib_node)
    {
        bool ok = true;

        for (pugi::xml_node node = lib_node.first_child(); node != nullptr && ok; node = node.next_sibling())
        {
            if(std::wstring(node.name()) == L"light")
            {
                LightSource *lgt;
                if((lgt = load_lightsource(node, scene.metadata)) == nullptr) return false;
                lgt->raw_xml = node;
                scene.light_sources[lgt->id] = lgt;
            }

        }

        return ok;
    }


    bool load_instanced_scene(InstancedScene &scene, pugi::xml_node scene_node)
    {
        bool ok = true;

        scene.custom_data = scene_node;
        scene.id = scene_node.attribute(L"id").as_uint(INVALID_ID);
        scene.name = ws2s(scene_node.attribute(L"name").as_string());
        if (!scene_node.attribute(L"bbox").empty())
        {
            scene.bbox = AABBFromString(scene_node.attribute(L"bbox").as_string());
        }

        //parse remap lists
        pugi::xml_node rl_nodes = scene_node.child(L"remap_lists");
        if (!rl_nodes.empty())
        {
            for (pugi::xml_node rl_node = rl_nodes.first_child(); rl_node != nullptr && ok; rl_node = rl_node.next_sibling())
            {
                uint32_t remap_list_id = rl_node.attribute(L"id").as_uint(INVALID_ID);
                if (remap_list_id == INVALID_ID)
                {
                    ok = false;
                    printf("[HydraScene::load_instanced_scene] Invalid remap list id\n");
                }
                else
                {
                    InstancedScene::RemapList remap_list;
                    remap_list.id = remap_list_id;
                    std::string rl_array_str = ws2s(rl_node.attribute(L"val").as_string());
                    auto indices = split(rl_array_str, ' ');
                    for (auto index : indices)
                    {
                        const char *str_c = index.c_str();
                        char *end_c = nullptr;
                        uint32_t value = strtol(index.c_str(), &end_c, 10);
                        if (end_c != str_c) //conversion happened
                            remap_list.remap.push_back((uint32_t)value);
                    }
                    if (remap_list.remap.size() == 0)
                    {
                        printf("[HydraScene::load_instanced_scene] Invalid remap list\n");
                        ok = false;
                    }
                    else
                        scene.remap_lists[remap_list_id] = remap_list;
                }
            }
        }

        if (!ok)
            return false;

        //parse instances and light instances
        for (pugi::xml_node inst_node = scene_node.first_child(); inst_node != nullptr && ok; inst_node = inst_node.next_sibling())
        {
            if (std::wstring(inst_node.name()) == L"instance")
            {
                Instance inst;
                inst.custom_data = inst_node;
                inst.id = inst_node.attribute(L"id").as_uint(INVALID_ID);
                inst.mesh_id = inst_node.attribute(L"mesh_id").as_uint(INVALID_ID);
                inst.rmap_id = inst_node.attribute(L"rmap_id").as_uint(INVALID_ID);
                inst.scn_id = inst_node.attribute(L"scn_id").as_uint(INVALID_ID);
                inst.scn_sid = inst_node.attribute(L"scn_sid").as_uint(INVALID_ID);
                inst.light_id = inst_node.attribute(L"light_id").as_uint(INVALID_ID);
                inst.linst_id = inst_node.attribute(L"linst_id").as_uint(INVALID_ID);

                if (inst.id == INVALID_ID || inst.mesh_id == INVALID_ID)
                {
                    printf("[HydraScene::load_instanced_scene] Invalid instance, each instance must have a unique id and a valid mesh (geom) id\n");
                    ok = false;
                }
                else if (inst_node.attribute(L"matrix").empty())
                {
                    printf("[HydraScene::load_instanced_scene] Invalid instance, each instance must have a transform matrix\n");
                    ok = false;
                }
                else
                {
                    inst.matrix = wstring_to_float4x4(inst_node.attribute(L"matrix").as_string());
                    scene.instances[inst.id] = inst;
                }
            }
            else if (std::wstring(inst_node.name()) == L"instance_light")
            {
                LightInstance inst;
                inst.custom_data = inst_node;
                inst.id = inst_node.attribute(L"id").as_uint(INVALID_ID);
                inst.mesh_id = inst_node.attribute(L"mesh_id").as_uint(INVALID_ID);
                inst.light_id = inst_node.attribute(L"light_id").as_uint(INVALID_ID);
                inst.lgroup_id = inst_node.attribute(L"lgroup_id").as_int(INVALID_ID); //it can be -1
                inst.matrix = wstring_to_float4x4(inst_node.attribute(L"matrix").as_string());

                if (inst.id == INVALID_ID || inst.light_id == INVALID_ID)
                {
                    printf("[HydraScene::load_instanced_scene] Invalid light instance, each light instance must have a unique id and a valid light id\n");
                    ok = false;
                }
                else if (inst_node.attribute(L"matrix").empty())
                {
                    printf("[HydraScene::load_instanced_scene] Invalid instance, each instance must have a transform matrix\n");
                    ok = false;
                }
                else
                {
                    inst.matrix = wstring_to_float4x4(inst_node.attribute(L"matrix").as_string());
                    scene.light_instances[inst.id] = inst;
                }
            }
        }

        if (scene.instances.size() == 0)
        {
            printf("[HydraScene::load_instanced_scene] No valid instances found on scene %d\n", scene.id);
            return false;
        }

        return true;
    }

    bool load_all_instanced_scenes(HydraScene &scene, pugi::xml_node lib_node)
    {
        bool ok = true;

        for (pugi::xml_node scene_node = lib_node.first_child(); scene_node != nullptr && ok; scene_node = scene_node.next_sibling())
        {
            if (std::wstring(scene_node.name()) == L"scene")
            {
                uint32_t id = scene_node.attribute(L"id").as_uint(INVALID_ID);
                if (id == INVALID_ID)
                {
                    printf("[HydraScene::load_instanced_scenes] Invalid scene id\n");
                    ok = false;
                }
                else
                {
                    InstancedScene inst_scene;
                    ok = load_instanced_scene(inst_scene, scene_node);
                    if (ok)
                        scene.scenes[id] = inst_scene;
                }
            }
        }

        return ok;
    }

    bool load_materials(HydraScene &scene, pugi::xml_node &lib_node);

    bool HydraScene::load(const std::string &filename)
    {
        clear();

        std::filesystem::path path(filename);
        if (!std::filesystem::exists(path))
        {
            printf("[HydraScene::load] Scene file %s does not exist\n", filename.c_str());
            return false;
        }
        metadata.scene_xml_path = filename;
        metadata.scene_xml_folder = path.parent_path().string();

        auto loaded = metadata.xml_doc.load_file(path.c_str());
        if (!loaded)
        {
            printf("[HydraScene::load] Failed to load xml file %s. xml_parse_status: %d\n", filename.c_str(), (int)loaded.status);
            return false;
        }

        pugi::xml_node root = metadata.xml_doc;
        if(metadata.xml_doc.child(L"root") != nullptr)
            root = metadata.xml_doc.child(L"root");

        metadata.custom_data = root;

        pugi::xml_node texturesLib  = root.child(L"textures_lib");
        pugi::xml_node materialsLib = root.child(L"materials_lib");
        pugi::xml_node geometryLib  = root.child(L"geometry_lib");
        pugi::xml_node lightsLib    = root.child(L"lights_lib");
        //pugi::xml_node spectraLib   = root.child(L"spectra_lib");

        pugi::xml_node cameraLib    = root.child(L"cam_lib");
        pugi::xml_node settingsNode = root.child(L"render_lib");
        pugi::xml_node scenesNode   = root.child(L"scenes");

        bool all_part_found = true;
        if (texturesLib.empty())
        {
            printf("[HydraScene::load] No textures lib\n");
            all_part_found = false;
        }
        if (materialsLib.empty())
        {
            printf("[HydraScene::load] No materials lib\n");
            all_part_found = false;
        }
        if (geometryLib.empty())
        {
            printf("[HydraScene::load] No geometry lib\n");
            all_part_found = false;
        }
        if (lightsLib.empty())
        {
            printf("[HydraScene::load] No lights lib\n");
            all_part_found = false;
        }
        // if (spectraLib.empty())
        // {
        //   printf("[HydraScene::load] No spectra lib\n");
        //   all_part_found = false;
        // }
        if (cameraLib.empty())
        {
            printf("[HydraScene::load] No camera lib\n");
            all_part_found = false;
        }
        if (settingsNode.empty())
        {
            printf("[HydraScene::load] No settings lib\n");
            all_part_found = false;
        }
        if (scenesNode.empty())
        {
            printf("[HydraScene::load] No scenes lib\n");
            all_part_found = false;
        }

        if (!all_part_found)
            return false;
        
        bool t_loaded = load_textures(*this, texturesLib);
        bool m_loaded = load_materials(*this, materialsLib);
        bool g_loaded = load_geometry(*this, geometryLib);
        bool l_loaded = load_lightsources(*this, lightsLib);
        //load_spectra(*this, spectraLib);
        bool c_loaded = load_cameras(*this, cameraLib);
        bool rs_loaded = load_all_render_settings(*this, settingsNode);
        bool s_loaded = load_all_instanced_scenes(*this, scenesNode);

        return t_loaded && m_loaded && g_loaded && l_loaded && c_loaded && rs_loaded && s_loaded;
    }

    bool save_geometry(const HydraScene &scene, const SceneMetadata &save_metadata, pugi::xml_node &lib_node)
    {
        for (const auto &[id, geom] : scene.geometries)
        {
            auto node = geom->custom_data ? lib_node.append_copy(geom->custom_data) : lib_node.append_child(L"geometry");
            geom->load_data(scene.metadata);
            geom->save_data(save_metadata);
            if(!geom->save_node(node)) return false;
        }
        return !lib_node.empty();
    }

    void set_texture(pugi::xml_node &colorNode, const TextureInstance &inst);
    void save_color_holder(pugi::xml_node &node, const ColorHolder &clr, bool allow_spectrum = true);

    bool save_lightsource(const LightSource *lgt, const SceneMetadata &newmeta, pugi::xml_node &node) {

        set_attr(node, L"id", lgt->id);
        set_attr(node, L"name", s2ws(lgt->name));
        set_attr(node, L"mat_id", lgt->mat_id);

        auto intNode = set_child(node, L"intensity");
        auto colNode = set_child(intNode, L"color");
        switch(lgt->type()) {
        case LightSource::Type::SKY:
            {
                set_attr(node, L"type", L"sky");
                //set_attr(node, L"shape", L"point"); is it necessary?
                set_attr(node, L"distribution", L"uniform");

                const LightSourceSky *sptr = static_cast<const LightSourceSky *>(lgt);
                if(sptr->texture) {
                    set_texture(colNode, *sptr->texture);
                }
                if(sptr->camera_back) {
                    auto backNode = set_child(node, L"back");
                    set_texture(backNode, *sptr->camera_back);
                }

            } break;
        case LightSource::Type::DIRECTIONAL:
            set_attr(node, L"type", L"directional"); 
            break;
        case LightSource::Type::RECT:
            {   
                set_attr(node, L"type", L"area"); 
                set_attr(node, L"shape", L"rect");
                auto sizeNode = set_child(node, L"size");
                set_attr(sizeNode, L"half_length", lgt->half_width);
                set_attr(sizeNode, L"half_length", lgt->half_length);
            }
            break;
        case LightSource::Type::DISK:
            set_attr(node, L"type", L"area");
            set_attr(node, L"shape", L"disk");
            set_attr(set_child(node, L"size"), L"radius", lgt->radius);
            break;
        case LightSource::Type::SPHERE:
            set_attr(node, L"type", L"area");
            set_attr(node, L"shape", L"sphere");
            set_attr(set_child(node, L"size"), L"radius", lgt->radius);
        case LightSource::Type::POINT:
            set_attr(node, L"type", L"point");
            set_attr(node, L"shape", L"point");
            if(lgt->distribution == LightSource::Dist::SPOT) {
                set_attr(node, L"distribution", L"spot");

                const LightSourceSpot *ptr = static_cast<const LightSourceSpot *>(lgt);
                set_attr(set_child(node, L"falloff_angle"), L"val", ptr->angle1);
                set_attr(set_child(node, L"falloff_angle2"), L"val", ptr->angle2);

                if(ptr->projective) {
                    auto &proj = *ptr->projective;
                    auto projNode = set_child(node, L"projective");
                    set_attr(set_child(projNode, L"fov"), L"val", proj.fov);
                    set_attr(set_child(projNode, L"nearClipPlane"), L"val", proj.nearClipPlane);
                    set_attr(set_child(projNode, L"farClipPlane"), L"val", proj.farClipPlane);

                    if(proj.texture) {
                        auto projTexNode = set_child(projNode, L"texture");
                        set_texture(projTexNode, *proj.texture);
                    }
                }
            }
            else {
                if(lgt->distribution == LightSource::Dist::LAMBERT) {
                    set_attr(node, L"distribution", "lambert");
                }
                else if(lgt->distribution == LightSource::Dist::OMNI) {
                    set_attr(node, L"distribution", "uniform");
                }
                else return false; //error
            }
            break;
        }

        save_color_holder(colNode, lgt->color, true);
        set_attr(set_child(intNode, L"multiplier"), L"val", lgt->power);

        if(lgt->ies) {
            auto iesNode = set_child(node, L"ies");
            const auto &ies = *lgt->ies; 

            fs::path path{ies.file_path}; //is always absolute
            fs::path abs_new_path = fs::path(newmeta.geometry_folder) / "ies" / path.filename();
            fs::path rel_new_path = get_relative_if_possible(fs::path(newmeta.scene_xml_folder), abs_new_path);
            fs::copy(path, abs_new_path, fs::copy_options::update_existing
                                                          | fs::copy_options::recursive);
            set_attr(iesNode, L"loc", s2ws(rel_new_path.string()));

            set_attr(iesNode, L"point_area", int(ies.point_area));
            if(ies.matrix) {
                set_attr(iesNode, L"matrix", LM_to_wstring(*ies.matrix));
            }
        }
        return true;
    }

    bool save_lightsources(const HydraScene &scene, const SceneMetadata &meta, pugi::xml_node &lib_node)
    {
        for (const auto &[id, lgt] : scene.light_sources)
        {
            pugi::xml_node node = lgt->raw_xml ? lib_node.append_copy(lgt->raw_xml) : lib_node.append_child(L"light");
            if(!save_lightsource(lgt, meta, node)) return false;
        }
        return !lib_node.empty();
    }

    bool save_instanced_scene(const InstancedScene &scene, pugi::xml_node &scene_node)
    {
        set_attr(scene_node, L"id", scene.id);
        std::wstring bbox_str = AABBToString(scene.bbox);
        set_attr(scene_node, L"bbox", bbox_str);
        set_attr(scene_node, L"name", s2ws(scene.name));

        if (scene.remap_lists.size() > 0)
        {
            pugi::xml_node all_remap_lists_node = scene_node.append_child(L"remap_lists");
            for (const auto &[id, remap_list] : scene.remap_lists)
            {
                std::wstring list_str;
                for (const auto &elem : remap_list.remap)
                    list_str += std::to_wstring(elem) + L" ";
                
                pugi::xml_node remap_list_node = all_remap_lists_node.append_child(L"remap_list");
                remap_list_node.append_attribute(L"id").set_value(id);
                remap_list_node.append_attribute(L"size").set_value(remap_list.remap.size());
                remap_list_node.append_attribute(L"val").set_value(list_str.c_str());
            }
        }

        if (scene.instances.size() > 0)
        {
            for (const auto &[id, instance] : scene.instances)
            {
                pugi::xml_node inst_node;
                if (instance.custom_data) {
                    inst_node = scene_node.append_copy(instance.custom_data);
                }
                else {
                    inst_node = scene_node.append_child(L"instance");
                }
                set_attr(inst_node, L"id", id);
                set_attr(inst_node, L"mesh_id", instance.mesh_id);
                set_attr(inst_node, L"matrix", float4x4ToString(instance.matrix));

                if (instance.rmap_id != INVALID_ID)
                    set_attr(inst_node, L"rmap_id", instance.rmap_id);
                if (instance.scn_id != INVALID_ID)
                    set_attr(inst_node, L"scn_id", instance.scn_id);
                if (instance.scn_sid != INVALID_ID)
                    set_attr(inst_node, L"scn_sid", instance.scn_sid);
                if (instance.light_id != INVALID_ID)
                    set_attr(inst_node, L"light_id", instance.light_id);
                if (instance.linst_id != INVALID_ID)
                    set_attr(inst_node, L"linst_id", instance.linst_id);
            }
        }

        if (scene.light_instances.size() > 0)
        {
            for (const auto &[id, linst] : scene.light_instances)
            {
                pugi::xml_node linst_node;
                if (linst.custom_data) {
                    linst_node = scene_node.append_copy(linst.custom_data);
                }
                else {
                    linst_node = scene_node.append_child(L"instance_light");
                }
                set_attr(linst_node, L"id", id);
                set_attr(linst_node, L"light_id", linst.light_id);
                set_attr(linst_node, L"matrix", float4x4ToString(linst.matrix));

                if (linst.mesh_id != INVALID_ID)
                    set_attr(linst_node, L"mesh_id", linst.mesh_id);
                
                if (linst.lgroup_id != INVALID_ID)
                    set_attr(linst_node, L"lgroup_id", linst.lgroup_id);
            }
        }

        return !scene_node.empty();
    }

    bool save_instanced_scenes(const HydraScene &scene, pugi::xml_node lib_node)
    {
        for (const auto &[id, inst_scene] : scene.scenes)
        {
            pugi::xml_node scene_node;
            if (inst_scene.custom_data) {
                scene_node = lib_node.append_copy(inst_scene.custom_data);
            }
            else {
                scene_node = lib_node.append_child(L"scene");
            }
            
            //clear all instances, light_instances and remap lists, as they will be saved below
            pugi::xml_node child_node = scene_node.first_child();
            while (!child_node.empty())
            {
                pugi::xml_node next_node = child_node.next_sibling();
                if (std::wstring(child_node.name()) == L"instance" || 
                        std::wstring(child_node.name()) == L"instance_light" ||
                        std::wstring(child_node.name()) == L"remap_lists")
                {
                    scene_node.remove_child(child_node);
                }
                child_node = next_node;
            }
            save_instanced_scene(inst_scene, scene_node);
        }
        return !lib_node.empty();
    }


    bool save_camera(const Camera &cam, pugi::xml_node &cam_node)
    {
        set_attr(cam_node, L"id", cam.id);
        set_attr(cam_node, L"name", s2ws(cam.name));

        set_child(cam_node, L"fov", cam.fov);
        set_child(cam_node, L"nearClipPlane", cam.nearPlane);
        set_child(cam_node, L"farClipPlane", cam.farPlane);


        if(cam.exposureMult != 1.0f) {
            set_child(cam_node, L"exposure_mult", cam.exposureMult);
        }

        set_child(cam_node, L"position", LM_to_wstring(cam.pos));
        set_child(cam_node, L"look_at", LM_to_wstring(cam.lookAt));
        set_child(cam_node, L"up", LM_to_wstring(cam.up));

        if(cam.has_matrix)
        {
            auto matrixNode = set_child(cam_node, L"matrix");
            set_attr(matrixNode, L"val", LM_to_wstring(LiteMath::transpose(cam.matrix)));
        }


        return true;
    }


    bool save_cameras(const HydraScene &scene, pugi::xml_node lib_node)
    {
        for (const auto &[id, cam] : scene.cameras)
        {
            pugi::xml_node cam_node;
            if (cam.custom_data) {
                cam_node = lib_node.append_copy(cam.custom_data);
            }
            else {
                cam_node = lib_node.append_child(L"camera");
            }
            if(!save_camera(cam, cam_node)) return false;
        }
        return !lib_node.empty();
    }


    bool save_render_settings(const RenderSettings &settings, pugi::xml_node &node)
    {
        set_attr(node, L"id", settings.id);
        set_attr(node, L"name", s2ws(settings.name));
        set_attr(node, L"type", L"HydraModern");


        set_child(node, L"width", std::to_wstring(settings.width));
        set_child(node, L"height", std::to_wstring(settings.height));
        set_child(node, L"trace_depth", std::to_wstring(settings.depth));
        set_child(node, L"diff_trace_depth", std::to_wstring(settings.depthDiffuse));
        set_child(node, L"maxRaysPerPixel", std::to_wstring(settings.spp));

        return true;
    }


    bool save_all_render_settings(const HydraScene &scene, pugi::xml_node lib_node)
    {
        for (const auto &[id, sett] : scene.render_settings)
        {
            pugi::xml_node sett_node;
            if (sett.custom_data) {
                sett_node = lib_node.append_copy(sett.custom_data);
            }
            else {
                sett_node = lib_node.append_child(L"render_settings");
            }
            if(!save_render_settings(sett, sett_node)) return false;
        }
        return !lib_node.empty();
    }


    bool save_textures(const HydraScene &scene, pugi::xml_node lib_node, const SceneMetadata &new_meta)
    {
        for (const auto &[id, tex] : scene.textures)
        {
            //auto tex_node = lib_node.append_copy(tex.custom_data);
            auto tex_node = lib_node.append_child(L"texture");
            if(!tex.save_info(tex_node, scene.metadata.scene_xml_folder, new_meta)) return false;
        }
        return !lib_node.empty();
    }

    bool save_materials(const HydraScene &scene, pugi::xml_node lib_node);

    bool HydraScene::save(const std::string &filename, const std::string &geometry_folder)
    {
        std::filesystem::file_status file_status = std::filesystem::status(filename);
        if (file_status.type() != std::filesystem::file_type::regular &&
                file_status.type() != std::filesystem::file_type::not_found)
        {
            printf("[HydraScene::save] File %s is a folder or some special file. It's name cannot be used as a scene file.\n", 
                         filename.c_str());
            return false;
        }
        std::filesystem::file_status status = std::filesystem::status(geometry_folder);
        if (status.type() == std::filesystem::file_type::not_found)
        {
            bool created = std::filesystem::create_directories(geometry_folder);
            if (!created)
            {
                printf("[HydraScene::save] Failed to create geometry folder %s\n", geometry_folder.c_str());
                return false;
            }
        }
        else if (status.type() != std::filesystem::file_type::directory)
        {
            printf("[HydraScene::save] Geometry folder %s is an existing file\n", geometry_folder.c_str());
            return false;
        }
        SceneMetadata save_metadata;
        save_metadata.scene_xml_path = filename;
        save_metadata.scene_xml_folder = std::filesystem::path(filename).parent_path().string();
        save_metadata.geometry_folder = geometry_folder;
        save_metadata.geometry_folder_relative = std::filesystem::relative(save_metadata.geometry_folder, save_metadata.scene_xml_folder).string();


        pugi::xml_document doc;
        pugi::xml_node texturesLib  = doc.append_child(L"textures_lib");
        pugi::xml_node materialsLib = doc.append_child(L"materials_lib");
        pugi::xml_node geometryLib  = doc.append_child(L"geometry_lib");
        pugi::xml_node lightsLib    = doc.append_child(L"lights_lib");
        //pugi::xml_node spectraLib   = root.child(L"spectra_lib");

        pugi::xml_node cameraLib    = doc.append_child(L"cam_lib");
        pugi::xml_node settingsLib = doc.append_child(L"render_lib");
        pugi::xml_node scenesLib   = doc.append_child(L"scenes");

         if (!save_textures(*this, texturesLib, save_metadata))
           return false;
        if (!save_materials(*this, materialsLib))
           return false;
        if (!save_geometry(*this, save_metadata, geometryLib))
            return false;
        if (!save_lightsources(*this, save_metadata, lightsLib))
           return false;
        if (!save_cameras(*this, cameraLib))
            return false;
        if (!save_all_render_settings(*this, settingsLib))
            return false;
        if (!save_instanced_scenes(*this, scenesLib))
            return false;

        metadata.geometry_folder = save_metadata.geometry_folder;
        metadata.geometry_folder_relative = save_metadata.geometry_folder_relative;
        metadata.scene_xml_folder = save_metadata.scene_xml_folder;
        metadata.scene_xml_path = save_metadata.scene_xml_path;

        return doc.save_file(filename.c_str());
    }

    unsigned HydraScene::get_total_number_of_primitives() const
    {
        unsigned count = 0;
        for (const auto &[id, geom] : geometries)
        {
            auto *mesh = dynamic_cast<const MeshGeometry *>(geom);
            if (mesh != nullptr) {
                if (mesh->is_loaded) {
                    count += mesh->mesh.TrianglesNum();
                }
                else {
                    count += mesh->custom_data.attribute(L"triNum").as_uint(0);
                }
            } 
            else {
                count += geom->custom_data.attribute(L"num_primitives").as_uint(0);
            }
        }
        return count;
    }

    uint32_t HydraScene::add_geometry(LiteScene::Geometry *geom)
    {
        uint32_t id = geometries.size();
        geom->id = id;
        geom->custom_data = metadata.xml_doc.child(L"geometry_lib").append_child(s2ws(geom->type_name).c_str());
        geometries[id] = geom;
        return id;
    }

    uint32_t HydraScene::add_mesh(const cmesh4::SimpleMesh &mesh)
    {
        uint32_t id = geometries.size();
        MeshGeometry *geom = new MeshGeometry();
        geom->mesh = mesh;
        geom->type_id = Geometry::MESH_TYPE_ID;
        geom->id = id;
        geom->name = "unnamed_mesh_"+std::to_string(id);
        geom->relative_file_path = "mesh_"+std::to_string(id) + ".vsgf";
        geom->type_name = "vsgf";
        geom->bytesize = mesh.SizeInBytes();
        geom->custom_data = metadata.xml_doc.child(L"geometry_lib").append_child(L"mesh");
        geom->custom_data.append_attribute(L"vertNum").set_value(mesh.VerticesNum());
        geom->custom_data.append_attribute(L"triNum").set_value(mesh.TrianglesNum());
        geom->is_loaded = true;
        geometries[id] = geom;
        return id;
    }
  
    uint32_t HydraScene::add_instance(uint32_t geomId, LiteMath::float4x4 transform)
    {
        auto it = geometries.find(geomId);
        if (it == geometries.end()) {
            printf("add_instance error - trying to add instance to geomId=%u that does not exist!\n", geomId);
            return uint32_t(-1);
        }
        else {
            if (scenes.empty()) {
                scenes[0] = InstancedScene();
                scenes[0].id = 0;
                scenes[0].bbox.boxMin = LiteMath::float3(-1,-1,-1);
                scenes[0].bbox.boxMax = LiteMath::float3( 1, 1, 1);
            }
            int instId = scenes[0].instances.size();
            scenes[0].instances[instId] = Instance();
            scenes[0].instances[instId].matrix = transform;
            scenes[0].instances[instId].mesh_id = geomId;
            return instId;
        }
    }


}