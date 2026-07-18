#include "integrator_pt_scene.h"
#include "neural_loader.h"
#include <cassert>
#include <utility>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

Material LoadNeuralBrdfMaterial(const std::string &scn_dir, const pugi::xml_node& materialNode,
                                std::vector<float> &m_neural_weights, std::vector<uint32_t> &m_neural_weights_offsets)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  uint32_t id = materialNode.attribute(L"id").as_uint();
  Material mat = {};
  mat.mtype = MAT_TYPE_NEURAL_BRDF;
  mat.lightId = uint(-1);

  const auto nnNode = materialNode.child(L"nn");

  //Loading weights
  std::string weights_path = fs::path(scn_dir) / hydra_xml::ws2s(nnNode.attribute(L"weights_loc").as_string());

  size_t weights_offset = m_neural_weights.size();

  nn::WeightsLoader wloader{weights_path};
  while(wloader.has_next())
  {
    const size_t next_size = wloader.next_size();
    const size_t old_size = m_neural_weights.size();
    m_neural_weights.resize(old_size + next_size);

    wloader.load_next(m_neural_weights.data() + old_size);
  }
  m_neural_weights_offsets[id] = weights_offset;
  return mat;
}

Material LoadKanBrdfMaterial(const std::string &scn_dir, const pugi::xml_node& materialNode,
                             std::vector<float> &m_neural_weights, std::vector<uint32_t> &m_neural_weights_offsets)
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

  //Loading latent texture
  /*std::vector<std::pair<HydraSampler, uint32_t>> loaded_tex = LoadLatentTexturesFromNode(texNode, texturesInfo, texCache, textures, scn_dir);
  m_neural_tex_offsets[id] = {m_neural_tex_ids.size(), loaded_tex.size()};
  for(const auto &[sampler_out, loaded_tex_id] : loaded_tex)
  {
    m_neural_tex_ids.push_back(loaded_tex_id);
  }*/


  //Loading weights
  std::string weights_path = fs::path(scn_dir) / hydra_xml::ws2s(nnNode.attribute(L"weights_loc").as_string());

  size_t weights_offset = m_neural_weights.size();

  nn::WeightsLoader wloader{weights_path};
  while(wloader.has_next())
  {
    const size_t next_size = wloader.next_size();
    const size_t old_size = m_neural_weights.size();
    m_neural_weights.resize(old_size + next_size);
    //std::vector<float> weights;
    //weights.resize(mat_size);
    //wloader.load_next(weights.data(), m_neural_weights.data() + old_size + mat_size);
    //nn::Transpose(weights.data(), m_neural_weights.data() + old_size, rows, cols);

    wloader.load_next(m_neural_weights.data() + old_size);
  }
  m_neural_weights_offsets[id] = weights_offset;
  return mat;
}