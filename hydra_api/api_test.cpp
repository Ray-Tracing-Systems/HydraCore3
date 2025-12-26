#include "hydra_api.h"
#include "mesh_utils.h"

#include <iostream>
#include <sstream>

int main(int argc, const char** argv)
{
  HR2_StorageRef scnStorage = hr2CreateStorage(HR2_STORAGE_CPU, HR2_ReserveOpions());

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // first load/create heavy scene data
  //

  HR2_CommandBuffer storageLevel = hr2CommandBufferStorage(scnStorage, HR2_CLEAR_AND_APPEND);

  // (1) Create materials
  //

  HR2_MaterialRef mat0 = hr2CreateMaterial(storageLevel);
  HR2_MaterialRef mat1 = hr2CreateMaterial(storageLevel);

  // set material #0
  {
    auto node = hr2MaterialParamNode(storageLevel, mat0);
    node.append_attribute(L"name") = L"MyFirstMaterial";
    node.append_attribute(L"type") = L"hydra_material";

    auto diff = node.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type") = L"lambert";
    
    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"0.5 0.75 0.5");
  }
  
  // set material #1
  {
    auto node = hr2MaterialParamNode(storageLevel, mat1);
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

  HR2_GeomRef cubeRef  = hr2CreateMeshFromData(storageLevel, "cube", inputCube);
  HR2_GeomRef planeRef = hr2CreateMeshFromData(storageLevel, "plane", inputPlane);
  
  // (3) Create lights
  //

  HR2_LightRef rectLight = hr2CreateLight(storageLevel);
  {
    auto lightNode = hr2LightParamNode(storageLevel, rectLight);
    
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

  HR2_SceneRef  sceneRef = hr2CreateScene(storageLevel);
  
  HR2_FrameBufferInfo fbInfo = {};
  fbInfo.width  = 512;
  fbInfo.height = 512;

  HR2_FrameImgRef frameImageRef = hr2CreateFrameImg(storageLevel, fbInfo);

  hr2Commit(storageLevel); // now scene library is finished 

  // MUST NOT USE storageLevel after hr2Commit, immediately report error! 
  // ALWAYS create new command buffer

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // next we can do relativelly quck scene create/update during edit/work with program
  //
  HR2_CommandBuffer sceneLvl = hr2CommandBufferScene(sceneRef, HR2_CLEAR_AND_APPEND); 

  // (4) Create camera; 
  //
  
  HR2_CameraRef camRef = hr2CreateCamera(sceneLvl);
  {
    auto camNode = hr2CameraParamNode(storageLevel, camRef);
    
    camNode.append_attribute(L"name")             = L"my_camera";
    camNode.append_child(L"fov").text()           = 45;
    camNode.append_child(L"nearClipPlane").text() = 0.01f;
    camNode.append_child(L"farClipPlane").text()  = 100.0;
    
    camNode.append_child(L"up").text().set(L"0 1 0");
    camNode.append_child(L"position").text().set(L"0 3 20");
    camNode.append_child(L"look_at").text().set(L"0 0 0");
  }
  
  // (5) render settings ... 
  //

  HR2_SettingsRef settingsRef = hr2CreateSettings(sceneLvl);
  {
    auto node = hr2SettingsParamNode(storageLevel, settingsRef);
    
    node.append_child(L"width").text()  = 512;
    node.append_child(L"height").text() = 512;
    
    node.append_child(L"method_primary").text()   = L"pathtracing";
    node.append_child(L"trace_depth").text()      = 6;
    node.append_child(L"diff_trace_depth").text() = 4;
    node.append_child(L"maxRaysPerPixel").text()  = 256;
    node.append_child(L"qmc_variant").text()      = 7; // enable all of them, results to '7'
  }

  hr2Commit(sceneLvl); // now scene is finished and we can render some scene

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // (6) Create scene as instances of existing objects and lights
  //
  const int NFrames = 10;
  for(int frame = 0; frame < NFrames; frame++)
  {
      
    // #NOTE: each frame we discard old scene and create the new one by instancing of existing geometry and lights
    //
    HR2_CommandBuffer frameLvl = hr2CommandBufferScene(sceneRef, HR2_CLEAR_AND_APPEND);
    {
    
      float m1[16] = {1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1,};
    
      float m2[16] = {1, 0, 0, 0,
                      0, 1, 0, 1,
                      0, 0, 1, 0,
                      0, 0, 0, 1,};
    
      float m3[16] = {1, 0, 0, 0,
                      0, 1, 0, 5,
                      0, 0, 1, 0,
                      0, 0, 0, 1,};
    
      hr2GeomInstance(frameLvl, planeRef, m1);
      hr2GeomInstance(frameLvl, cubeRef, m2);
      hr2LightInstance(frameLvl, rectLight, m3);
      
      for (int z = -2; z <= 2; z++)
      {
        for (int x = -2; x <= 2; x++)
        {
          float m4[16] = {1, 0, 0, float(x) * float(frame),
                          0, 1, 0, 1,
                          0, 0, 1, float(z) * float(frame),
                          0, 0, 0, 1,};
        
          hr2GeomInstance(frameLvl, cubeRef, m4);
        }
      }
    }

    hr2CommitAndRender(frameLvl, camRef, settingsRef, frameImageRef, false); 
    
    // when async commit-and-wait is used
    //
    //while (true)
    //{
    //  std::this_thread::sleep_for(std::chrono::milliseconds(100));
    //  auto info = hr2HaveUpdate(frameLvl);
    //
    //  if (info.haveUpdateFB)
    //  {
    //    std::cout << "rendering progress = " << info.progress << "% \r";
    //    std::cout.flush();
    //  }
    //
    //  if (info.finalUpdate)
    //    break;
    //}

    // render to frameImageRef

    std::stringstream stream;
    stream << "z_img_" << frame << ".png";
    std::string fileName = stream.str();
   
    std::cout << "save frame " << frame << " to " << fileName.c_str() << std::endl;
    hr2SaveFrameBuffer(frameImageRef, fileName.c_str());
    
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  return 0;
}