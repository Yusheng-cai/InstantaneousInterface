/* This file is part of PyMesh. Copyright (c) 2015 by Qingnan Zhou */
#pragma once

#include <vector>
#include "tools/Algorithm.h"

#include "Mesh.h"
#include "tools/CommonTypes.h"

class LongEdgeRemoval {
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        LongEdgeRemoval(Mesh& mesh) :
            m_vertices_(mesh.getvertices()), m_faces_(mesh.gettriangles()), mesh_(mesh) {}

    public:
        void run(Real max_length, bool recursive=true);

        std::vector<vertex> get_vertices() const { return m_vertices_; }
        std::vector<triangle> get_faces() const { return m_faces_; }

        VectorI get_ori_faces() const { return m_ori_faces; }

    private:
        void init_edge_map();
        void split_long_edges(Real max_length);
        size_t retriangulate();
        void triangulate_chain(
                std::vector<VectorI>& faces,
                const std::vector<size_t>& chain,
                size_t v0_idx, size_t v1_idx, size_t v2_idx);
        std::vector<size_t> get_vertex_chain_around_triangle(
                size_t fi, size_t& v0_idx, size_t& v1_idx, size_t& v2_idx);

    private:
        Mesh& mesh_;

        std::vector<vertex> m_vertices_;
        std::vector<triangle> m_faces_;
        VectorI  m_ori_faces;

        DupletMap<size_t> m_edge_map;
};