#include "hydra_api.h"
#include "mesh_utils.h"

int main(int argc, const char** argv)
{
  HAPI_SceneLibrary  scnLibrary   = hapiCreateLibraryEmpty(HAPI_STORAGE_CPU, HAPI_ReserveOpions());
  HAPI_CommandBuffer appendBuffer = hapiCreateCommandBuffer(scnLibrary, HAPI_CMD_APPEND);

  HAPI_Material mat0 = hapiCreateMaterialEmpty(appendBuffer, "MyFirstMaterial");
  HAPI_Material mat1 = hapiCreateMaterialEmpty(appendBuffer, "MySecondMaterial");

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
  
  SimpleMesh meshPlane = CreatePlane(10.0f);
  SimpleMesh meshCube  = CreateCube(1.0f);
  
  HAPI_Geom cubeRef = hapiCreateMeshEmpty(appendBuffer, "cube");
  {
    //hrMeshVertexAttribPointer4f(cubeRef, L"pos",      &meshCube.vPos[0]);
    //hrMeshVertexAttribPointer4f(cubeRef, L"norm",     &meshCube.vNorm[0]);
    //hrMeshVertexAttribPointer2f(cubeRef, L"texcoord", &meshCube.vTexCoord[0]);
    //
    //hrMeshMaterialId(cubeRef, mat0.id);
    //
    //hrMeshAppendTriangles3(cubeRef, int(meshCube.triIndices.size()), &meshCube.triIndices[0]);
  }
  
  HAPI_Geom planeRef = hapiCreateMeshEmpty(appendBuffer, "plane");
  {
    //hrMeshVertexAttribPointer4f(planeRef, L"pos",      &meshPlane.vPos[0]);
    //hrMeshVertexAttribPointer4f(planeRef, L"norm",     &meshPlane.vNorm[0]);
    //hrMeshVertexAttribPointer2f(planeRef, L"texcoord", &meshPlane.vTexCoord[0]);
    //
    //hrMeshMaterialId(planeRef, mat1.id);
    //
    //hrMeshAppendTriangles3(planeRef, int(meshPlane.triIndices.size()), &meshPlane.triIndices[0]);
  }
 
  hapiCommitCommandBuffer(appendBuffer); // now scene library is finished and we can render some scene
  

  return 0;
}