#include "integrator_pt_scene.h"
#include <iostream>
#include <filesystem>

static void LoadMerlMaterial(const std::string &path, Material &mat,
                             std::vector<float> &a_measured_brdfs, std::vector<Integrator::MeasuredBrdfEntry> &a_measured_brdf_data)
{

}

Material LoadMeasuredMaterial(const std::string &scn_dir,
                              const pugi::xml_node& materialNode,
                              std::vector<float> &a_measured_brdfs, std::vector<Integrator::MeasuredBrdfEntry> &a_measured_brdf_data)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  uint32_t id = materialNode.attribute(L"id").as_uint();
  Material mat = {};
  mat.mtype = MAT_TYPE_MEASURED;
  mat.lightId = uint(-1);

  const auto dataNode = materialNode.child(L"data");
  
  
  
  const std::wstring type = dataNode.attribute(L"type").as_string();
  const std::filesystem::path data_path = std::filesystem::path(scn_dir) / hydra_xml::ws2s(dataNode.attribute(L"loc").as_string());

  if(type == L"merl" || data_path.extension() == ".binary") {
    LoadMerlMaterial(data_path, mat, a_measured_brdfs, a_measured_brdf_data);
  }
  else {
    std::cout << "[LoadMeasuredMaterial] Unknown measured material type: " + hydra_xml::ws2s(type) << std::endl;

    Integrator::MeasuredBrdfEntry entry;
    entry.dim = {0, 0, 0, 0};
    entry.offset = 0;
    a_measured_brdf_data.push_back(std::move(entry));
  }


  return mat;
}