#include "include/cmaterial.h"
#include "integrator_pt_scene.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <vector>

static bool safe_read_exact(std::istream &file, char *dst, uint64_t bytecount, std::string &errmesg)
{
  try {
    file.read(dst, bytecount);
  }
  catch(std::ios::failure &e) {
    errmesg = e.code().message();
    return false;
  }

  if(file.gcount() != bytecount) {
    errmesg = "Unexpected EOF";
    return false;
  }

  return true;
}

static void LoadMerlMaterial(const std::string &path, Material &mat,
                             std::vector<float> &a_measured_brdfs, std::vector<Integrator::MeasuredBrdfEntry> &a_measured_brdf_data)
{
  constexpr int SAMPLING_THETA_H = 90;
  constexpr int SAMPLING_THETA_D = 90;
  constexpr int SAMPLING_PHI_D = 180; //360 is halved in merl code
  constexpr uint64_t MERL_SIZE = SAMPLING_PHI_D * SAMPLING_THETA_D * SAMPLING_THETA_H;
  constexpr double RED_SCALE = 1.0;
  constexpr double GREEN_SCALE = 1.15;
  constexpr double BLUE_SCALE = 1.66;
  constexpr double SCALE_DIV = 1500.0;

  std::string errmesg;
  std::ifstream file{path, std::ios::in | std::ios::binary};

  try {
    file.exceptions(std::ios::badbit | std::ios::failbit);
  }
  catch(std::ios::failure &e) {
    errmesg = e.code().message();
    throw std::runtime_error("Error opening file (" + path + "): " + errmesg);
  }

  uint32_t dims[3];
  if(!safe_read_exact(file, reinterpret_cast<char *>(dims), 3 * sizeof(uint32_t), errmesg)) {
     throw std::runtime_error("Error reading MERL BRDF file (" + path + "): " + errmesg);
  }

  if(dims[0] != SAMPLING_THETA_H && dims[1] != SAMPLING_THETA_D && dims[2] != SAMPLING_PHI_D) {
    throw std::runtime_error("Incorrect data shape in MERL BRDF file (" + path + ")");
  }

  
  std::vector<double> data64(MERL_SIZE * 3);

  if(!safe_read_exact(file, reinterpret_cast<char *>(data64.data()), MERL_SIZE * sizeof(double) * 3, errmesg)) {
     throw std::runtime_error("Error reading MERL BRDF file (" + path + "): " + errmesg);
  }
  file.close();


  Integrator::MeasuredBrdfEntry entry;
  entry.dim = uint4(SAMPLING_THETA_H, SAMPLING_THETA_D, 1, SAMPLING_PHI_D);
  entry.offset = uint32_t(a_measured_brdfs.size());
  entry.nchannels = 3;
  entry.phi_range = M_PI;

  a_measured_brdfs.resize(entry.offset + MERL_SIZE * 3);
  float *data32 = a_measured_brdfs.data() + entry.offset;

  for(uint64_t i = 0; i < MERL_SIZE; ++i) {
    data32[i * 3 + 0] = max(0.0f, static_cast<float>(data64[i]                 * RED_SCALE / SCALE_DIV));
    data32[i * 3 + 1] = max(0.0f, static_cast<float>(data64[i + MERL_SIZE]     * GREEN_SCALE / SCALE_DIV));
    data32[i * 3 + 2] = max(0.0f, static_cast<float>(data64[i + 2 * MERL_SIZE] * BLUE_SCALE / SCALE_DIV));

    //std::cout << data32[i * 3 + 0] << " " << data32[i * 3 + 1] << " " << data32[i * 3 + 2] << std::endl;
  }

  mat.datai[MEASURED_DATAIDX] = static_cast<uint>(a_measured_brdf_data.size());
  a_measured_brdf_data.push_back(std::move(entry));

}

static void LoadHydraMeasuredMaterial(const std::string &path, Material &mat,
                                      std::vector<float> &a_measured_brdfs, std::vector<Integrator::MeasuredBrdfEntry> &a_measured_brdf_data)
{
  static constexpr char HEADER_STRING[] = "hydrameasured1";

  std::string errmesg;
  std::ifstream file{path, std::ios::in | std::ios::binary};

  try {
    file.exceptions(std::ios::badbit | std::ios::failbit);
  }
  catch(std::ios::failure &e) {
    errmesg = e.code().message();
    throw std::runtime_error("Error opening file (" + path + "): " + errmesg);
  }
  
  char buf[sizeof(HEADER_STRING)];
  file.read(buf, sizeof(HEADER_STRING) - 1);
  buf[sizeof(HEADER_STRING) - 1] = '\0';
  if(strcmp(buf, HEADER_STRING) != 0) {
    throw std::runtime_error("File does not contain nn weigths: " + path);
  }

  uint32_t dims[4];
  if(!safe_read_exact(file, reinterpret_cast<char *>(dims), 4 * sizeof(uint32_t), errmesg)) {
     throw std::runtime_error("Error reading Hydra Measured BRDF file (" + path + "): " + errmesg);
  }
  uint8_t n_channels;
  if(!safe_read_exact(file, reinterpret_cast<char *>(&n_channels), 1, errmesg)) {
     throw std::runtime_error("Error reading Hydra Measured BRDF file (" + path + "): " + errmesg);
  }
  float phi_range;
  if(!safe_read_exact(file, reinterpret_cast<char *>(&phi_range), sizeof(float), errmesg)) {
     throw std::runtime_error("Error reading Hydra Measured BRDF file (" + path + "): " + errmesg);
  }


  uint64_t data_size = dims[0] * dims[1] * dims[2] * dims[3];
  //std::cout << dims[0] << " " << dims[1] << " " << dims[2] << " " << dims[3] << " " << uint(n_channels) << " " << std::endl;

  Integrator::MeasuredBrdfEntry entry;
  entry.dim = uint4(dims[0], dims[1], dims[2], dims[3]);
  entry.offset = uint32_t(a_measured_brdfs.size());
  entry.nchannels = n_channels;
  entry.phi_range = phi_range;

  a_measured_brdfs.resize(entry.offset + data_size * n_channels);
  float *data32 = a_measured_brdfs.data() + entry.offset;

  if(!safe_read_exact(file, reinterpret_cast<char *>(data32), data_size * n_channels * sizeof(float), errmesg)) {
     throw std::runtime_error("Error reading Hydra Measured BRDF file (" + path + "): " + errmesg);
  }
  file.close();

  mat.datai[MEASURED_DATAIDX] = static_cast<uint>(a_measured_brdf_data.size());
  a_measured_brdf_data.push_back(std::move(entry));
}

Material LoadMeasuredMaterial(const std::string &scn_dir,
                              const pugi::xml_node& materialNode,
                              std::vector<float> &a_measured_brdfs, std::vector<Integrator::MeasuredBrdfEntry> &a_measured_brdf_data)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  //uint32_t id = materialNode.attribute(L"id").as_uint();
  Material mat = {};
  mat.mtype = MAT_TYPE_MEASURED;
  mat.lightId = uint(-1);

  const auto dataNode = materialNode.child(L"data");
  const float alpha = hydra_xml::readval1f(materialNode.child(L"alpha"), 0.1f);
  mat.data[MEASURED_ALPHA] = alpha;
  
  
  const std::wstring type = dataNode.attribute(L"type").as_string(L"");
  const std::filesystem::path data_path = std::filesystem::path(scn_dir) / hydra_xml::ws2s(dataNode.attribute(L"loc").as_string());

  if(type == L"merl" || (type.empty() && data_path.extension() == ".binary")) {
    LoadMerlMaterial(data_path, mat, a_measured_brdfs, a_measured_brdf_data);
  }
  else if(type == L"hydra" || (type.empty() && data_path.extension() == ".hydram")) {
    LoadHydraMeasuredMaterial(data_path, mat, a_measured_brdfs, a_measured_brdf_data);
  }
  else {
    std::cout << "[LoadMeasuredMaterial] Unknown measured material type: " + hydra_xml::ws2s(type) << std::endl;
    mat.datai[MEASURED_DATAIDX] = uint(-1);
  }


  return mat;
}