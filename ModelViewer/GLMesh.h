#pragma once

#include "tools/CommonTypes.h"
#include "src/Mesh.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glad.h"

#include <vector>
#include <array>

namespace ModelView 
{
    struct GLVertex {
        // position
        glm::vec3 Position;
        glm::vec3 Normal;
    };

    class GLMesh{
        public:
            GLMesh(std::string name);
            void setupMesh();
        

        private:
            // input filename --> ply file etc.
            std::string filename_;
            std::vector<GLVertex> vertices_;
            std::vector<unsigned int> indices_;

            int nf_, nv_;

            Mesh m_;

            // buffer data 
            unsigned int VBO_, EBO_, VAO_;
    };

};