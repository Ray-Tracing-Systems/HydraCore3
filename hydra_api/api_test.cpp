#include "hydra_api.h"
#include "mesh_utils.h"

int main(int argc, const char** argv)
{
  HR2_SceneLibraryRef scnLibrary = hr2CreateLibrary(HR2_STORAGE_CPU, HR2_ReserveOpions());
  HR2_CommandBuffer appendBuffer = hr2CreateCommandBuffer(scnLibrary, HR2_CMDBUF_APPEND);

  // (1) Create materials
  //

  HR2_MaterialRef mat0 = hr2CreateMaterial(appendBuffer);
  HR2_MaterialRef mat1 = hr2CreateMaterial(appendBuffer);

  // set material #0
  {
    auto node = hr2MaterialParamNode(mat0);
    node.append_attribute(L"name") = L"MyFirstMaterial";
    node.append_attribute(L"type") = L"hydra_material";

    auto diff = node.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type") = L"lambert";
    
    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"0.5 0.75 0.5");
  }
  
  // set material #1
  {
    auto node = hr2MaterialParamNode(mat1);
    node.append_attribute(L"name") = L"MySecondMaterial";
    node.append_attribute(L"type") = L"hydra_material";
    
    auto diff = node.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type") = L"lambert";
    
    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"0.25 0.25 0.25");
  }
  
  // (1) Create geometry
  //

  SimpleMesh meshPlane = CreatePlane(10.0f);
  SimpleMesh meshCube  = CreateCube(1.0f);
  
  HR2_MeshInput inputCube, inputPlane;
  {
    inputCube.vPosPtr       = meshCube.vPos.data();
    inputCube.vNormPtr      = meshCube.vNorm.data();
    inputCube.vTexCoordPtr  = meshCube.vTexCoord.data();
    inputCube.matIdAll      = mat0.id;
    inputCube.indicesPtr    = (uint32_t*)meshCube.triIndices.data();
    inputCube.indicesNum    = uint32_t(meshCube.triIndices.size());
  
    inputPlane.vPosPtr      = meshPlane.vPos.data();
    inputPlane.vNormPtr     = meshPlane.vNorm.data();
    inputPlane.vTexCoordPtr = meshPlane.vTexCoord.data();
    inputPlane.matIdAll     = mat1.id;
    inputPlane.indicesPtr   = (uint32_t*)meshPlane.triIndices.data();
    inputPlane.indicesNum   = uint32_t(meshPlane.triIndices.size());
  }

  HR2_GeomRef cubeRef  = hr2CreateMeshFromData(appendBuffer, "cube", inputCube);
  HR2_GeomRef planeRef = hr2CreateMeshFromData(appendBuffer, "plane", inputPlane);
  
  // (3) Create lights
  //

  HR2_LightRef rectLight = hr2CreateLight(appendBuffer);
  {
    auto lightNode = hr2LightParamNode(rectLight);
    
    lightNode.append_attribute(L"name").set_value(L"my_area_light");
    lightNode.append_attribute(L"type").set_value(L"area");
    lightNode.append_attribute(L"shape").set_value(L"rect");
    lightNode.append_attribute(L"distribution").set_value(L"diffuse"); // you can use both 'set_value' or '='
    lightNode.append_attribute(L"visible").set_value(L"1");
    
    auto sizeNode = lightNode.append_child(L"size");
    
    sizeNode.append_attribute(L"half_length") = 1.0f;
    sizeNode.append_attribute(L"half_width")  = 1.0f;
    
    auto intensityNode = lightNode.append_child(L"intensity");
    
    intensityNode.append_child(L"color").append_attribute(L"val")      = L"1 1 1";
    intensityNode.append_child(L"multiplier").append_attribute(L"val") = 25.0f;
  }

  // (4) Create camera; TODO: do we need to create camera via command buffer ?
  //
  
  HR2_CameraRef camRef = hr2CreateCamera(appendBuffer);
  {
    auto camNode = hr2CameraParamNode(camRef);
    
    camNode.append_attribute(L"name")             = L"my_camera";
    camNode.append_child(L"fov").text()           = 45;
    camNode.append_child(L"nearClipPlane").text() = 0.01f;
    camNode.append_child(L"farClipPlane").text()  = 100.0;
    
    camNode.append_child(L"up").text().set(L"0 1 0");
    camNode.append_child(L"position").text().set(L"0 3 20");
    camNode.append_child(L"look_at").text().set(L"0 0 0");
  }
 
  hr2CommitCommandBuffer(appendBuffer); // now scene library is finished and we can render some scene
  
  // (5) render settings ... TODO: do we need to create camera via command buffer ?
  //

  HR2_SettingsRef settingsRef = hr2CreateSettings(appendBuffer);
  {
    auto node = hr2SettingsParamNode(settingsRef);
    
    node.append_child(L"width").text()  = 512;
    node.append_child(L"height").text() = 512;
    
    node.append_child(L"method_primary").text()   = L"pathtracing";
    node.append_child(L"trace_depth").text()      = 6;
    node.append_child(L"diff_trace_depth").text() = 4;
    node.append_child(L"maxRaysPerPixel").text()  = 256;
    node.append_child(L"qmc_variant").text()      = 7; // enable all of them, results to '7'
  }

  // (6) Create scene as instances of existing objects and lights
  //

  return 0;
}