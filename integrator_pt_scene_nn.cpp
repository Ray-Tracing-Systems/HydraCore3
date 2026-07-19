#include "integrator_pt_scene.h"
#include "neural_loader.h"
#include <cassert>
#include <utility>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

static void LoadNeuralWeights(uint32_t mat_id, const std::string &scn_dir, const pugi::xml_node& nnNode,
                              std::vector<float> &neural_weights, std::vector<uint32_t> &neural_weights_offsets)
{
  std::string weights_path = fs::path(scn_dir) / hydra_xml::ws2s(nnNode.attribute(L"weights_loc").as_string());

  size_t weights_offset = neural_weights.size();

  nn::WeightsLoader wloader{weights_path};
  while(wloader.has_next())
  {
    const size_t next_size = wloader.next_size();
    const size_t old_size = neural_weights.size();
    neural_weights.resize(old_size + next_size);

    wloader.load_next(neural_weights.data() + old_size);
  }
  neural_weights_offsets[mat_id] = uint32_t(weights_offset);
}

Material LoadNeuralBrdfMaterial(const std::string &scn_dir, const pugi::xml_node& materialNode,
                                std::vector<float> &neural_weights, std::vector<uint32_t> &neural_weights_offsets)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  uint32_t id = materialNode.attribute(L"id").as_uint();
  Material mat = {};
  mat.mtype = MAT_TYPE_NEURAL_BRDF;
  mat.lightId = uint(-1);

  const auto medianGridNode = materialNode.child(L"median");
  uint32_t median_id = medianGridNode.attribute(L"mat_id").as_uint();
  mat.datai[NBRDF_MEDIANIDX] = median_id;

  const auto nnNode = materialNode.child(L"nn");

  LoadNeuralWeights(id, scn_dir, nnNode, neural_weights, neural_weights_offsets);
  return mat;
}

Material LoadKanBrdfMaterial(const std::string &scn_dir, const pugi::xml_node& materialNode,
                             std::vector<float> &neural_weights, std::vector<uint32_t> &neural_weights_offsets)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  uint32_t id = materialNode.attribute(L"id").as_uint();
  Material mat = {};
  mat.mtype = MAT_TYPE_KANBRDF;
  mat.lightId = uint(-1);
  //TODO
  const float alpha = hydra_xml::readval1f(materialNode.child(L"alpha"), 0.1f);
  mat.data[KANBRDF_ALPHA] = alpha;

  const auto nnNode = materialNode.child(L"nn");

  LoadNeuralWeights(id, scn_dir, nnNode, neural_weights, neural_weights_offsets);
  
  return mat;
}