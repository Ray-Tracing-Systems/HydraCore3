#include <LiteScene/scene.h>
#include "hydraxml.h"

namespace {

    HydraSampler ReadSamplerFromColorNode(const pugi::xml_node a_colorNodes, bool from_spectrum)
    {
        hydra_xml::HydraSampler res;
        pugi::xml_node texNode;

        if(from_spectrum)
            texNode = a_colorNodes.child(L"spectrum");
        else
            texNode = a_colorNodes.child(L"texture");
        if(texNode == nullptr)
            return res;
        
        res.texId = texNode.attribute(L"id").as_uint();
        
        if(texNode.attribute(L"addressing_mode_u") != nullptr)
        {
            std::wstring addModeU = texNode.attribute(L"addressing_mode_u").as_string();
            res.sampler.addressU  = GetAddrModeFromString(addModeU);
        } 

        if(texNode.attribute(L"addressing_mode_v") != nullptr)
        {
            std::wstring addModeV = texNode.attribute(L"addressing_mode_v").as_string();
            res.sampler.addressV  = GetAddrModeFromString(addModeV);
        }

        if(texNode.attribute(L"addressing_mode_w") == nullptr)
            res.sampler.addressW  = res.sampler.addressV;
        else
        {
            std::wstring addModeW = texNode.attribute(L"addressing_mode_w").as_string();
            res.sampler.addressW  = GetAddrModeFromString(addModeW);
        }

        res.sampler.filter = Sampler::Filter::LINEAR;
        if(texNode.attribute(L"filter") != nullptr)
        {
            std::wstring filterMode = texNode.attribute(L"filter").as_string();
            if(filterMode == L"point" || filterMode == L"nearest")
                res.sampler.filter = Sampler::Filter::NEAREST;
            else if(filterMode == L"cubic" || filterMode == L"bicubic")
                res.sampler.filter = Sampler::Filter::CUBIC;
        }

        if(texNode.attribute(L"input_gamma") != nullptr)
            res.inputGamma = texNode.attribute(L"input_gamma").as_float();

        const std::wstring inputAlphaMode = texNode.attribute(L"input_alpha").as_string();
        if(inputAlphaMode == L"alpha")
            res.alphaFromRGB = false;
        
        // read texture matrix
        //
        std::wstringstream inputStream(texNode.attribute(L"matrix").as_string()); // in HydraXML we store matrices by rows
        for(int i=0;i<4;i++)
            inputStream >> res.row0[i];
        for(int i=0;i<4;i++)
            inputStream >> res.row1[i];
        return res;
    }

}

namespace ls {



}