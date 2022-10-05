#include "GLMesh.h"

ModelView::GLMesh::GLMesh(std::string filename)
:filename_(filename)
{
    // read mesh 
    MeshTools::readPLYlibr(filename_, m_);

    // copy data
    const auto& f = m_.gettriangles();
    const auto& v  = m_.getvertices();
    nf_ = f.size();
    nv_ = v.size(); 


    // copy triangle data
    indices_.resize(nf_*3);
    for (int i=0;i<nf_;i++){
        for (int j=0;j<3;j++){
            indices_[i*3+j] = f[i][j];
        }
    }

    // copy vertex data 
    vertices_.resize(nv_);
    for (int i=0;i<nv_;i++){
        glm::vec3 vec(v[i][0], v[i][1], v[i][2]);
        glm::vec3 n(v[i].normals_[0], v[i].normals_[1], v[i].normals_[2]);
        GLVertex glV = {vec,n};
        vertices_[i] = glV;
    }
}

void ModelView::GLMesh::setupMesh(){
    glGenVertexArrays(1, &VAO_);
    glGenBuffers(1, &VBO_);
    glGenBuffers(1, &EBO_);

    glBindVertexArray(VAO_);
}