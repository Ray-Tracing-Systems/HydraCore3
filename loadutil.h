#ifndef LITESCENE_LOADUTIL_H_
#define LITESCENE_LOADUTIL_H_
#include "3rd_party/pugixml.hpp"
#include "hydraxml.h"
#include "LiteMath.h"
#include <type_traits>
#include <locale>
#include <codecvt>
#include <vector>
#include <sstream>
#include <filesystem>
#include <numeric>

namespace fs = std::filesystem;

namespace LiteScene {

    inline pugi::xml_node set_child(pugi::xml_node node, const pugi::char_t *name)
    {
        pugi::xml_node child = node.child(name) ? node.child(name) : node.append_child(name);
        return child;
    }

    template<typename T>
    inline pugi::xml_node set_child(pugi::xml_node node, const pugi::char_t *name, T value)
    {
        auto child = set_child(node, name);
        child.text().set(value);
        return child;
    }

    inline pugi::xml_node set_child(pugi::xml_node node, const pugi::char_t *name, const std::wstring &value)
    {
        return set_child<const pugi::char_t *>(node, name, value.c_str());
    }


    template<typename T>
    inline void set_attr(pugi::xml_node node, const pugi::char_t *name, T value)
    {
        if (node.attribute(name).empty()) node.append_attribute(name).set_value(value); 
        else node.attribute(name).set_value(value); 
    }

    inline void set_attr(pugi::xml_node node, const pugi::char_t *name, const std::wstring &value)
    {
        set_attr<const pugi::char_t *>(node, name, value.c_str());
    }


    inline std::wstring LM_to_wstring(const LiteMath::float3 &v)
    {
        return std::to_wstring(v.x) + L" " + std::to_wstring(v.y) + L" " + std::to_wstring(v.z);
    }

    inline std::wstring LM_to_wstring(const LiteMath::float4 &v)
    {
        return std::to_wstring(v.x) + L" " + std::to_wstring(v.y) + L" " + std::to_wstring(v.z) + L" " + std::to_wstring(v.w);
    }

    inline std::wstring LM_to_wstring(const LiteMath::float4x4 &v)
    {
        return LM_to_wstring(v.get_row(0)) + L" "
             + LM_to_wstring(v.get_row(1)) + L" " 
             + LM_to_wstring(v.get_row(2)) + L" " 
             + LM_to_wstring(v.get_row(3));
    }
/*
    inline LiteMath::float3 to_float3(const std::wstring &str)
    {
        LiteMath::float3 res;
        std::wstringstream ss{str};
        ss >> res.x;
    }
*/

    inline std::vector<float> wstring_to_float_arr(const std::wstring &str, int count)
    {
        std::vector<float> result(count);
        std::wstringstream inputStream(str);
        for (int i = 0; i < count; i++)
        {
            inputStream >> result[i];
        }
        return result;
    }

    inline LiteMath::float4x4 wstring_to_float4x4(const std::wstring &str)
    {
        auto data = wstring_to_float_arr(str, 16);
        LiteMath::float4x4 result;

        result.set_row(0, LiteMath::float4(data[0],data[1], data[2], data[3]));
        result.set_row(1, LiteMath::float4(data[4],data[5], data[6], data[7]));
        result.set_row(2, LiteMath::float4(data[8],data[9], data[10], data[11]));
        result.set_row(3, LiteMath::float4(data[12],data[13], data[14], data[15])); 

        return result;
    }

    inline float pop_attr_float(pugi::xml_node &node, const pugi::char_t *name)
    {
        auto val = node.attribute(name).as_float();
        node.remove_attribute(name);
        return val;
    }

    inline uint32_t pop_attr_uint(pugi::xml_node &node, const pugi::char_t *name)
    {
        auto val = node.attribute(name).as_uint();
        node.remove_attribute(name);
        return val;
    }

    inline std::wstring pop_attr_str(pugi::xml_node &node, const pugi::char_t *name)
    {
        auto val = node.attribute(name).as_string();
        node.remove_attribute(name);
        return std::wstring(val);
    }


    //https://stackoverflow.com/questions/60530422/how-to-check-a-file-is-contained-in-a-folder-with-c
    // =======
    inline fs::path normalized_trimed(const fs::path& p)
    {
        auto r = p.lexically_normal();
        if (r.has_filename()) return r;
        return r.parent_path();
    }

    inline fs::path get_relative_if_possible(const fs::path& dir, const fs::path& target)
    {   
        if(!dir.empty()) {
            auto base = normalized_trimed(dir);
            auto sub = normalized_trimed(target).parent_path();
            auto m = std::mismatch(base.begin(), base.end(), 
                                   sub.begin(), sub.end());

            if(m.first == base.end()) { // dir contains target
                fs::path rel_path = std::accumulate(m.second, sub.end(), fs::path{}, std::divides{}) / target.filename();
                return rel_path.lexically_normal();
            }
        }
        return target.lexically_normal();
    }
    // =======
}
#endif