#include "integrator_pt_scene.h"
#include "neural.h"
#include <algorithm>
#include <cassert>
#include <utility>

//???
std::vector<std::pair<HydraSampler, uint32_t>> LoadLatentTexturesFromNode(const pugi::xml_node texNode, 
                                const std::vector<TextureInfo> &texturesInfo,
                                std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache,
                                std::vector<std::shared_ptr<ICombinedImageSampler>> &textures)
{

  const auto refs_attr = texNode.attribute(L"ids");
  std::vector<uint32_t> tex_refs = hydra_xml::readvalVectorU(refs_attr);

  Sampler::AddressMode modeU, modeV, modeW;
  modeU = modeV = modeW = Sampler::AddressMode::WRAP;
  Sampler::Filter filter = Sampler::Filter::LINEAR;
  float4 row0{1,0,0,0};
  float4 row1{0,1,0,0};

  if(texNode.attribute(L"addressing_mode_u") != nullptr)
  {
    std::wstring addModeU = texNode.attribute(L"addressing_mode_u").as_string();
    modeU  = GetAddrModeFromString(addModeU);
  } 

  if(texNode.attribute(L"addressing_mode_v") != nullptr)
  {
    std::wstring addModeV = texNode.attribute(L"addressing_mode_v").as_string();
    modeW  = GetAddrModeFromString(addModeV);
  }

  if(texNode.attribute(L"addressing_mode_w") == nullptr)
    modeW  = modeV;
  else
  {
    std::wstring addModeW = texNode.attribute(L"addressing_mode_w").as_string();
    modeW  = GetAddrModeFromString(addModeW);
  }

  if(texNode.attribute(L"filter") != nullptr)
  {
    std::wstring filterMode = texNode.attribute(L"filter").as_string();
    if(filterMode == L"point" || filterMode == L"nearest")
      filter = Sampler::Filter::NEAREST;
    else if(filterMode == L"cubic" || filterMode == L"bicubic")
      filter = Sampler::Filter::CUBIC;
  }

  std::wstringstream inputStream(texNode.attribute(L"matrix").as_string()); // in HydraXML we store matrices by rows
  for(int i=0;i<4;i++)
    inputStream >> row0[i];
  for(int i=0;i<4;i++)
    inputStream >> row1[i];


  std::vector<std::pair<HydraSampler, uint32_t>> res;
  res.reserve(tex_refs.size());
  for(uint32_t texId : tex_refs)
  {
    HydraSampler s;
    s.sampler.addressU = modeU;
    s.sampler.addressV = modeV;
    s.sampler.addressW = modeW;

    s.sampler.filter = filter;
    s.row0 = row0;
    s.row1 = row1;
    s.texId = texId;
    const auto& [sampler_out, loaded_tex_id] = LoadTextureById(texId, texturesInfo, s, texCache, textures);
    res.push_back({s, loaded_tex_id});
  }

  return res;
}


Material LoadNeuralBrdfMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache,
                                std::vector<std::shared_ptr<ICombinedImageSampler>> &textures,
                                std::vector<uint> &m_neural_tex_ids, std::vector<uint2> &m_neural_tex_offsets,
                                std::vector<float> &m_neural_weights, std::vector<uint> &m_neural_weights_offsets)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  uint32_t id = materialNode.attribute(L"id").as_uint();
  Material mat = {};
  mat.mtype = MAT_TYPE_NEURAL_BRDF;
  mat.lightId = uint(-1);
  //TODO


  const auto nnNode = materialNode.child(L"nn");
  const auto texNode = nnNode.child(L"texture");

  //Loading latent texture
  std::vector<std::pair<HydraSampler, uint32_t>> loaded_tex = LoadLatentTexturesFromNode(texNode, texturesInfo, texCache, textures);
  m_neural_tex_offsets[id] = {m_neural_tex_ids.size(), loaded_tex.size()};
  for(const auto &[sampler_out, loaded_tex_id] : loaded_tex)
  {
    m_neural_tex_ids.push_back(loaded_tex_id);
  }


  //Loading weights
  std::string weights_path = hydra_xml::ws2s(nnNode.attribute(L"weights_loc").as_string());

  size_t weights_offset = m_neural_weights.size();
  nn::WeightsLoader wloader{weights_path};
  while(wloader.has_next())
  {
    size_t mat_size = wloader.next_rows() * wloader.next_cols();
    size_t old_size = m_neural_weights.size();
    m_neural_weights.resize(old_size + mat_size + wloader.next_rows());
    wloader.load_next(m_neural_weights.data() + old_size, m_neural_weights.data() + old_size + mat_size);
  }
  m_neural_weights_offsets[id] = weights_offset;


  return mat;
}