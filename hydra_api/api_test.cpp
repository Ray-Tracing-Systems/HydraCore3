#include "hydra_api.h"
#include "mesh_utils.h"

int main(int argc, const char** argv)
{
  HAPI_SceneLibrary  scnLibrary   = hapiCreateLibraryEmpty(HAPI_STORAGE_CPU, HAPI_ReserveOpions());
  HAPI_CommandBuffer appendBuffer = hapiCreateCommandBuffer(scnLibrary, HAPI_CMDBUF_APPEND);

  // (1) Create materials
  //

  HAPI_Material mat0 = hapiCreateMaterial(appendBuffer, "MyFirstMaterial");
  HAPI_Material mat1 = hapiCreateMaterial(appendBuffer, "MySecondMaterial");

  // set material #0
  {
    auto node = hapiMaterialParamNode(mat0);
    node.append_attribute(L"type") = L"???";

    auto diff = node.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type") = L"lambert";
    
    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"0.5 0.75 0.5");
  }
  
  // set material #1
  {
    auto node = hapiMaterialParamNode(mat1);
    node.append_attribute(L"type") = L"???";
    
    auto diff = node.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type") = L"lambert";
    
    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"0.25 0.25 0.25");
  }
  
  // (1) Create geometry
  //

  SimpleMesh meshPlane = CreatePlane(10.0f);
  SimpleMesh meshCube  = CreateCube(1.0f);
  
  HAPI_MeshInput inputCube, inputPlane;
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

  HAPI_Geom cubeRef  = hapiCreateMeshFromData(appendBuffer, "cube", inputCube);
  HAPI_Geom planeRef = hapiCreateMeshFromData(appendBuffer, "plane", inputPlane);
  
  // (3) Create lights
  //

  HAPI_Light rectLight = hapiCreateLight(appendBuffer, "my_area_light");
  {
    auto lightNode = hapiLightParamNode(rectLight);
    
    lightNode.append_attribute(L"type").set_value(L"area");
    lightNode.append_attribute(L"shape").set_value(L"rect");
    lightNode.append_attribute(L"distribution").set_value(L"diffuse"); // you can use both 'set_value' or '='
    
    auto sizeNode = lightNode.append_child(L"size");
    
    sizeNode.append_attribute(L"half_length") = 1.0f;
    sizeNode.append_attribute(L"half_width")  = 1.0f;
    
    auto intensityNode = lightNode.append_child(L"intensity");
    
    intensityNode.append_child(L"color").append_attribute(L"val")      = L"1 1 1";
    intensityNode.append_child(L"multiplier").append_attribute(L"val") = 25.0f;
  }

  // (4) Create camera; TODO: do we need to create camera via command buffer ?
  //
  
  HAPI_Camera camRef = hapiCreateCamera(appendBuffer, "my_camera");
  {
    auto camNode = hapiCameraParamNode(camRef);
    
    camNode.append_child(L"fov").text().set(L"45");
    camNode.append_child(L"nearClipPlane").text().set(L"0.01");
    camNode.append_child(L"farClipPlane").text().set(L"100.0");
    
    camNode.append_child(L"up").text().set(L"0 1 0");
    camNode.append_child(L"position").text().set(L"0 3 20");
    camNode.append_child(L"look_at").text().set(L"0 0 0");
  }
 
  hapiCommitCommandBuffer(appendBuffer); // now scene library is finished and we can render some scene
  
  // (5) render settings ... ?
  //


  // (6) Create scene as instances of existing objects and lights
  //

  return 0;
}