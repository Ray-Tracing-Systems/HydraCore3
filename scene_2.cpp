#include "scene_2.h"
#include "hydraxml.h"

#include <sstream>
#include <fstream>
#include <locale>
#include <codecvt>
#include <filesystem>

namespace LiteScene
{
  #define SET_ATTR(node, name, value) \
    {\
    if (node.attribute(name).empty()) node.append_attribute(name).set_value(value); \
    else node.attribute(name).set_value(value); \
    }

  pugi::xml_node set_child(pugi::xml_node &node, const pugi::char_t *name, const std::wstring &value)
  {
    pugi::xml_node child = node.child(name).empty() ? node.child(name) : node.append_child(name);
    child.text().set(value.c_str());
    return child;
  }

  pugi::xml_node set_child(pugi::xml_node &node, const pugi::char_t *name)
  {
    pugi::xml_node child = node.child(name).empty() ? node.child(name) : node.append_child(name);
    return child;
  }

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

  std::wstring LM_to_wstring(const LiteMath::float3 &v)
  {
    return std::to_wstring(v.x) + L" " + std::to_wstring(v.y) + L" " + std::to_wstring(v.z);
  }

  std::wstring LM_to_wstring(const LiteMath::float4 &v)
  {
    return std::to_wstring(v.x) + L" " + std::to_wstring(v.y) + L" " + std::to_wstring(v.z) + L" " + std::to_wstring(v.w);
  }

  std::wstring LM_to_wstring(const LiteMath::float4x4 &v)
  {
    return LM_to_wstring(v.get_row(0)) + L" "
         + LM_to_wstring(v.get_row(1)) + L" " 
         + LM_to_wstring(v.get_row(2)) + L" " 
         + LM_to_wstring(v.get_row(3));
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

  std::vector<float> load_float_array_from_string(const std::wstring &str, int count)
  {
    std::vector<float> result(count);
    std::wstringstream inputStream(str);
    for (int i = 0; i < count; i++)
    {
      inputStream >> result[i];
    }
    return result;
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

  LiteMath::float4x4 float4x4FromString(const std::wstring &matrix_str)
  {
    auto data = load_float_array_from_string(matrix_str, 16);
    LiteMath::float4x4 result;

    result.set_row(0, LiteMath::float4(data[0],data[1], data[2], data[3]));
    result.set_row(1, LiteMath::float4(data[4],data[5], data[6], data[7]));
    result.set_row(2, LiteMath::float4(data[8],data[9], data[10], data[11]));
    result.set_row(3, LiteMath::float4(data[12],data[13], data[14], data[15])); 

    return result;
  }

  std::wstring float4x4ToString(const LiteMath::float4x4 &matrix)
  {
    return save_float_array_to_string({ matrix.get_row(0).x, matrix.get_row(0).y, matrix.get_row(0).z, matrix.get_row(0).w,
                                       matrix.get_row(1).x, matrix.get_row(1).y, matrix.get_row(1).z, matrix.get_row(1).w,
                                       matrix.get_row(2).x, matrix.get_row(2).y, matrix.get_row(2).z, matrix.get_row(2).w,
                                       matrix.get_row(3).x, matrix.get_row(3).y, matrix.get_row(3).z, matrix.get_row(3).w });
  }

  AABB AABBFromString(const std::wstring &str)
  {
    auto data = load_float_array_from_string(str, 6);
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
  
  pugi::xml_node Geometry::save_node_base() const
  {
    pugi::xml_node node = custom_data;
    
    SET_ATTR(node, L"id", id);
    SET_ATTR(node, L"bytesize", bytesize);
    SET_ATTR(node, L"name", s2ws(name).c_str());
    SET_ATTR(node, L"type", s2ws(type_name).c_str());

    return node;
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

  pugi::xml_node MeshGeometry::save_node() const
  {
    if (relative_file_path == INVALID_PATH)
    {
      printf("[MeshGeometry::save_node] No location is specified. Save data first\n");
      return pugi::xml_node();
    }
    pugi::xml_node node = save_node_base();
    node.set_name(L"mesh");
    SET_ATTR(node, L"loc", s2ws(relative_file_path).c_str());
    return node;
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

  pugi::xml_node CustomGeometry::save_node() const
  {
    return save_node_base();
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

  void HydraScene::clear()
  {
    for (int i = 0; i < textures.size(); i++)
      delete textures[i];

    for (int i = 0; i < materials.size(); i++)
      delete materials[i];

    for (int i = 0; i < geometries.size(); i++)
      delete geometries[i];

    metadata = SceneMetadata();

    textures.clear();
    materials.clear();
    geometries.clear();

    light_sources.clear();
    cameras.clear();
    render_settings.clear();
    scenes.clear();
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
        cam.matrix = LiteMath::transpose(float4x4FromString(cam_node.child(L"matrix").attribute(L"val").as_string()));
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
          inst.matrix = float4x4FromString(inst_node.attribute(L"matrix").as_string());
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
        inst.matrix = float4x4FromString(inst_node.attribute(L"matrix").as_string());

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
          inst.matrix = float4x4FromString(inst_node.attribute(L"matrix").as_string());
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
    
    //load_textures(*this, texturesLib);
    //load_materials(*this, materialsLib);
    bool g_loaded = load_geometry(*this, geometryLib);
    //load_lights(*this, lightsLib);
    //load_spectra(*this, spectraLib);
    bool c_loaded = load_cameras(*this, cameraLib);
    bool rs_loaded = load_all_render_settings(*this, settingsNode);
    bool s_loaded = load_all_instanced_scenes(*this, scenesNode);

    return g_loaded && c_loaded && rs_loaded && s_loaded;
  }

  bool save_geometry(const HydraScene &scene, const SceneMetadata &save_metadata, pugi::xml_node lib_node)
  {
    for (const auto &[id, geom] : scene.geometries)
    {
      geom->load_data(scene.metadata);
      geom->save_data(save_metadata);
      pugi::xml_node geom_node = geom->save_node();
      lib_node.append_copy(geom_node);
    }
    return !lib_node.empty();
  }

  bool save_instanced_scene(const InstancedScene &scene, pugi::xml_node &scene_node)
  {
    SET_ATTR(scene_node, L"id", scene.id);
    std::wstring bbox_str = AABBToString(scene.bbox);
    SET_ATTR(scene_node, L"bbox", bbox_str.c_str());
    SET_ATTR(scene_node, L"name", s2ws(scene.name).c_str());

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
        pugi::xml_node inst_node = scene_node.append_copy(instance.custom_data);
        SET_ATTR(inst_node, L"id", id);
        SET_ATTR(inst_node, L"mesh_id", instance.mesh_id);
        SET_ATTR(inst_node, L"matrix", float4x4ToString(instance.matrix).c_str());

        if (instance.rmap_id != INVALID_ID)
          SET_ATTR(inst_node, L"rmap_id", instance.rmap_id);
        if (instance.scn_id != INVALID_ID)
          SET_ATTR(inst_node, L"scn_id", instance.scn_id);
        if (instance.scn_sid != INVALID_ID)
          SET_ATTR(inst_node, L"scn_sid", instance.scn_sid);
        if (instance.light_id != INVALID_ID)
          SET_ATTR(inst_node, L"light_id", instance.light_id);
        if (instance.linst_id != INVALID_ID)
          SET_ATTR(inst_node, L"linst_id", instance.linst_id);
      }
    }

    if (scene.light_instances.size() > 0)
    {
      for (const auto &[id, linst] : scene.light_instances)
      {
        pugi::xml_node linst_node = scene_node.append_copy(linst.custom_data);
        SET_ATTR(linst_node, L"id", id);
        SET_ATTR(linst_node, L"light_id", linst.light_id);
        SET_ATTR(linst_node, L"matrix", float4x4ToString(linst.matrix).c_str());

        if (linst.mesh_id != INVALID_ID)
          SET_ATTR(linst_node, L"mesh_id", linst.mesh_id);
        
        if (linst.lgroup_id != INVALID_ID)
          SET_ATTR(linst_node, L"lgroup_id", linst.lgroup_id);
      }
    }

    return !scene_node.empty();
  }

  bool save_instanced_scenes(const HydraScene &scene, pugi::xml_node lib_node)
  {
    for (const auto &[id, inst_scene] : scene.scenes)
    {
      auto scene_node = lib_node.append_copy(inst_scene.custom_data);
      
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
      SET_ATTR(cam_node, L"id", cam.id);
      SET_ATTR(cam_node, L"name", s2ws(cam.name).c_str());

      set_child(cam_node, L"fov", std::to_wstring(cam.fov));
      set_child(cam_node, L"nearClipPlane", std::to_wstring(cam.nearPlane));
      set_child(cam_node, L"farClipPlane", std::to_wstring(cam.farPlane));


      if(cam.exposureMult != 1.0f) {
        set_child(cam_node, L"exposure_mult", std::to_wstring(cam.exposureMult));
      }

      set_child(cam_node, L"position", LM_to_wstring(cam.pos));
      set_child(cam_node, L"look_at", LM_to_wstring(cam.lookAt));
      set_child(cam_node, L"up", LM_to_wstring(cam.up));

      if(cam.has_matrix)
      {
        auto matrixNode = set_child(cam_node, L"matrix");
        SET_ATTR(matrixNode, L"val", LM_to_wstring(LiteMath::transpose(cam.matrix)).c_str());
      }


      return true;
  }


  bool save_cameras(const HydraScene &scene, pugi::xml_node lib_node)
  {
    for (const auto &[id, cam] : scene.cameras)
    {
      auto cam_node = lib_node.append_copy(cam.custom_data);
      if(!save_camera(cam, cam_node)) return false;
    }
    return !lib_node.empty();
  }


  bool save_render_settings(const RenderSettings &settings, pugi::xml_node &node)
  {
    SET_ATTR(node, L"id", settings.id);
    SET_ATTR(node, L"name", s2ws(settings.name).c_str());
    SET_ATTR(node, L"type", L"HydraModern");


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
      auto sett_node = lib_node.append_copy(sett.custom_data);
      if(!save_render_settings(sett, sett_node)) return false;
    }
    return !lib_node.empty();
  }


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
      bool created = std::filesystem::create_directory(geometry_folder);
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

    // if (!save_textures(*this, texturesLib))
    //   return false;
    // if (!save_materials(*this, materialsLib))
    //   return false;
    if (!save_geometry(*this, save_metadata, geometryLib))
      return false;
    // if (!save_lights(*this, lightsLib))
    //   return false;
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


}