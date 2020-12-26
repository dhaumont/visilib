/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2017 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Visilib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Visilib. If not, see <http://www.gnu.org/licenses/>
*/

#pragma once

#include <vector>
#include <stdlib.h>
#include <unordered_map>

#include "math_vector_3.h"
#include "helper_triangle_mesh.h"
#include "helper_triangle_mesh_container.h"
#include "math_matrix_4.h"
#include "math_vector_2.h"

namespace visilib
{
    class HelperTriangleMesh;

    /** @brief An helper class allowing to generate procedural triangle mesh such as grid, cube, sphere, etc.
    
    This class is mainly used when writing test or demo code, to generate different occluders
    */
    class HelperSyntheticMeshBuilder
    {
    public:
        inline static size_t getKey(size_t myFaceNumber, size_t v0, size_t v1)
        {
            return v0 < v1 ? v0 * myFaceNumber + v1 : v1 * myFaceNumber + v0;
        }
        enum SubdivisionPatternType
        {
            Regular,
            OneInside,
            ThreeOutside
        };
    

        /** @brief Generate a regular grid of edge size 1 in the xy plane centered around the origin*/
        static HelperTriangleMesh* generateRegularGrid(int subdivision);

        /** @brief Generate a unit cube centered around the origin*/
        static HelperTriangleMesh* generateCube(int subdivision);

        /** @brief Generate a unit sphere centered around the origin*/
        static HelperTriangleMesh* generateSphere(int subdivision);

        /** @brief Generate a checkboard grid of edge size 1 in the xy plane centered around the origin*/
        static HelperTriangleMesh* generateCheckBoard(int subdivision, bool even);

        /** @brief Generate a slot in the xy plane centered around the origin*/
        static HelperTriangleMesh* generateSlot(float x0, float y0, float xSize, float ySize);

        /** @brief Generate the vertices of a regular polygon*/
        static void generateRegularPolygon(std::vector<MathVector3f>& vertices, size_t vertexCount);
   
        /** @brief Generate the vertices of a slot*/
        static void generateSlot(std::vector<int>& indices, std::vector<MathVector3f>& vertices, float x0, float y0,float xSize, float ySize);

        /** @brief Generate the vertex table and index table of a unit cube*/
        static void generateCube(std::vector<int>& indices, std::vector<MathVector3f>& vertices);

        /** @brief Removes in a random way a given number of faces*/
        static void removeFaces(HelperTriangleMesh* mesh, size_t count);

        /** @brief Add a random perturbation to a given mesh*/
        static void addRandomness(HelperTriangleMesh* mesh, float randomness);

        /** @brief Scales the input scene*/
        static void scale(HelperTriangleMeshContainer* scene, float scale);

        /** @brief Rotates the input scene*/
        static void rotate(HelperTriangleMeshContainer* scene, float pitch, float yaw, float roll);

        /** @brief Translated the input scene*/
        static void translate(HelperTriangleMeshContainer* scene, const MathVector3f& translation);

        /** @brief  Scales the input mesh*/
        static void scale(HelperTriangleMesh* mesh, float scale);

        /** @brief Rotates the input mesh*/
        static void rotate(HelperTriangleMesh* mesh, float pitch, float yaw, float roll);

        /** @brief Translates the input mesh*/
        static void translate(HelperTriangleMesh* mesh, const MathVector3f& translation);

        /** @brief Add a random perturbation to a given vertex table*/
        static void addRandomness(std::vector<MathVector3f>& vertices, float randomness);

        /** @brief  Scales the input vertex table*/
        static void scale(std::vector<MathVector3f>& vertices, float scale);

        /** @brief  Rotates the input vertex table*/
        static void rotate(std::vector<MathVector3f>& vertices, float pitch, float yaw, float roll);

        /** @brief  Translates the input vertex table*/
        static void translate(std::vector<MathVector3f>& vertices, const MathVector3f& translation);

        /** @brief Scale the input scene such that it fits in a unit cube*/
        static void rescaleToUnitBox(HelperTriangleMeshContainer* scene);
    private:

        /** @brief Subdivide a mesh given by its vertex table and index table following a given pattern*/
        static void subdivide(std::vector<int>& indices, std::vector<MathVector3f>& vertices, SubdivisionPatternType pattern = Regular);
    };

    inline void HelperSyntheticMeshBuilder::subdivide(std::vector<int>& indices, std::vector<MathVector3f>& vertices, SubdivisionPatternType pattern)
    {
        size_t size = indices.size();

        std::vector<std::pair<size_t, size_t>> s;
        s.push_back(std::pair<size_t, size_t>(0, 1));
        s.push_back(std::pair<size_t, size_t>(1, 2));
        s.push_back(std::pair<size_t, size_t>(0, 2));

        std::unordered_map<size_t, size_t> myEdgesMap;

        size_t i[6];
        indices.resize(size * 4);
        size_t lastIndex = size * 4 - 1;

        for (int index = (int)size - 3; index >= 0; index -= 3)
        {
            i[0] = indices[index + 0];
            i[1] = indices[index + 1];
            i[2] = indices[index + 2];
            i[3] = 0;
            i[4] = 0;
            i[5] = 0;

            for (size_t split = 0; split < 3; split++)
            {
                size_t a = s[split].first;
                size_t b = s[split].second;
                size_t key = getKey(size, i[a], i[b]);
                auto iter = myEdgesMap.find(key);
                if (iter != myEdgesMap.end())
                {
                    i[split + 3] = iter->second;
                }
                else
                {
                    i[split + 3] = vertices.size();
                    MathVector3f newVertex = (vertices[i[a]] + vertices[i[b]]) * 0.5;
                    vertices.push_back(newVertex);
                    myEdgesMap.insert(std::pair<size_t, size_t>(key, i[split + 3]));
                }
            }

            if (pattern == Regular)
            {
                indices[lastIndex--] = (int)i[0]; indices[lastIndex--] = (int)i[3];	indices[lastIndex--] = (int)i[5];
                indices[lastIndex--] = (int)i[3]; indices[lastIndex--] = (int)i[1];	indices[lastIndex--] = (int)i[4];
                indices[lastIndex--] = (int)i[5]; indices[lastIndex--] = (int)i[4];	indices[lastIndex--] = (int)i[2];
                indices[lastIndex--] = (int)i[3]; indices[lastIndex--] = (int)i[4];	indices[lastIndex--] = (int)i[5];
            }
            /*	else if (pattern == OneInside)
            {
            indices[index + 0] = i[3]; indices[index + 1] = i[4];	indices[index + 2] = i[5];
            }
            else if (pattern == ThreeOutside)
            {
            indices[index + 0] = i[0]; indices[index + 1] = i[3];	indices[index + 2] = i[5];
            indices.push_back(i[3]);   indices.push_back(i[1]);     indices.push_back(i[4]);
            indices.push_back(i[5]);   indices.push_back(i[4]);     indices.push_back(i[2]);
            }*/
        }
    }

    inline void HelperSyntheticMeshBuilder::generateCube(std::vector<int> & indices, std::vector<MathVector3f> & vertices)
    {
        float s = 0.5f;

        vertices.push_back(MathVector3f(-s, -s, -s));
        vertices.push_back(MathVector3f(-s, s, -s));
        vertices.push_back(MathVector3f(s, s, -s));
        vertices.push_back(MathVector3f(s, -s, -s));
        vertices.push_back(MathVector3f(-s, -s, s));
        vertices.push_back(MathVector3f(-s, s, s));
        vertices.push_back(MathVector3f(s, s, s));
        vertices.push_back(MathVector3f(s, -s, s));

        indices.push_back(0);   indices.push_back(1);   indices.push_back(2);
        indices.push_back(0);   indices.push_back(2);   indices.push_back(3);
        indices.push_back(4);   indices.push_back(6);   indices.push_back(5);
        indices.push_back(4);   indices.push_back(7);   indices.push_back(6);
        indices.push_back(4);   indices.push_back(5);   indices.push_back(1);
        indices.push_back(4);   indices.push_back(1);   indices.push_back(0);
        indices.push_back(3);   indices.push_back(2);   indices.push_back(6);
        indices.push_back(3);   indices.push_back(6);   indices.push_back(7);
        indices.push_back(1);   indices.push_back(5);   indices.push_back(6);
        indices.push_back(1);   indices.push_back(6);   indices.push_back(2);
        indices.push_back(4);   indices.push_back(0);   indices.push_back(3);
        indices.push_back(4);   indices.push_back(3);   indices.push_back(7);
    }

    inline HelperTriangleMesh * HelperSyntheticMeshBuilder::generateRegularGrid(int subdivision)
    {
        std::vector<MathVector3f> vertices;
        std::vector<int> indices;

        float s = 0.5f;

        vertices.push_back(MathVector3f(-s, -s, 0));
        vertices.push_back(MathVector3f(-s, s, 0));
        vertices.push_back(MathVector3f(s, s, 0));
        vertices.push_back(MathVector3f(s, -s, 0));

        indices.push_back(0); indices.push_back(1); indices.push_back(2);
        indices.push_back(0); indices.push_back(2); indices.push_back(3);

        for (size_t i = 0; i < subdivision; i++)
        {
            subdivide(indices, vertices);
        }

        return new HelperTriangleMesh(vertices, indices);
    }

    inline HelperTriangleMesh* HelperSyntheticMeshBuilder::generateCheckBoard(int subdivision, bool even)
    {
        auto result = generateRegularGrid(subdivision - 1);

        subdivide(result->getIndices(), result->getVertices(), even ? OneInside : ThreeOutside);

        return result;
    }

    inline HelperTriangleMesh * HelperSyntheticMeshBuilder::generateCube(int subdivision)
    {
        std::vector<MathVector3f> vertices;
        std::vector<int> indices;

        generateCube(indices, vertices);

        for (size_t i = 0; i < subdivision; i++)
        {
            subdivide(indices, vertices);
        }

        return new HelperTriangleMesh(vertices, indices);
    }

    inline HelperTriangleMesh* HelperSyntheticMeshBuilder::generateSphere(int subdivision)
    {
        std::vector<MathVector3f> vertices;
        std::vector<int> indices;

        generateCube(indices, vertices);

        for (size_t i = 0; i < subdivision; i++)
        {
            subdivide(indices, vertices);
        }

        for (size_t i = 0; i < vertices.size(); i++)
        {
            vertices[i].normalize();
        }

        return new HelperTriangleMesh(vertices, indices);
    }

    inline HelperTriangleMesh* HelperSyntheticMeshBuilder::generateSlot(float x0, float y0, float xSize, float ySize)
    {
        std::vector<MathVector3f> vertices;
        std::vector<int> indices;

        generateSlot(indices, vertices, x0, y0,xSize, ySize);

        return new HelperTriangleMesh(vertices, indices);
    }

    inline void HelperSyntheticMeshBuilder::generateSlot(std::vector<int>& indices, std::vector<MathVector3f>& vertices, float x0, float y0, float xSize, float ySize)
    {
        float s = 0.5f;
        float x = 0.5f * xSize;
        float y = 0.5f * ySize;

        vertices.push_back(MathVector3f(-s,  s, 0));
        vertices.push_back(MathVector3f( 0,  s, 0));
        vertices.push_back(MathVector3f( s,  s, 0));
        vertices.push_back(MathVector3f( s,  0, 0));
        vertices.push_back(MathVector3f( s, -s, 0));
        vertices.push_back(MathVector3f( 0, -s, 0));
        vertices.push_back(MathVector3f(-s, -s, 0));
        vertices.push_back(MathVector3f(-s,  0, 0));

        vertices.push_back(MathVector3f(x0-x, y0+y, 0));
        vertices.push_back(MathVector3f(x0+x, y0+y, 0));
        vertices.push_back(MathVector3f(x0+x, y0-y, 0));
        vertices.push_back(MathVector3f(x0-x, y0-y, 0));
     
        indices.push_back(0); indices.push_back(1);     indices.push_back(8);
        indices.push_back(1); indices.push_back(9);     indices.push_back(8);
        indices.push_back(1); indices.push_back(2);     indices.push_back(9);
        indices.push_back(9); indices.push_back(2);     indices.push_back(3);
        indices.push_back(9); indices.push_back(3);     indices.push_back(10);
        indices.push_back(10); indices.push_back(3);     indices.push_back(4);
        indices.push_back(10); indices.push_back(4);     indices.push_back(5);
        indices.push_back(11); indices.push_back(10);     indices.push_back(5);
        indices.push_back(6); indices.push_back(11);     indices.push_back(5);
        indices.push_back(7); indices.push_back(11);     indices.push_back(6);
        indices.push_back(7); indices.push_back(8);     indices.push_back(11);
        indices.push_back(7); indices.push_back(0);     indices.push_back(8);


    }

    inline void HelperSyntheticMeshBuilder::generateRegularPolygon(std::vector<MathVector3f> & vertices, size_t vertexCount)
    {
        if (vertexCount == 1)
        {
            vertices.push_back(MathVector3f::Zero());
        }
        else if (vertexCount == 2)
        {
            vertices.push_back(MathVector3f(-1, 0, 0));
            vertices.push_back(MathVector3f(1, 0, 0));
        }
        else
        {
            const float PI = 3.14159265358979323846f;

            float theta = 0;
            float dt = -2 * PI / (float)vertexCount;

            for (size_t i = 0; i < vertexCount; i++)
            {
                vertices.push_back(MathVector3f(cos(theta), sin(theta), 0));

                theta += dt;
            }
        }
    }

    inline void HelperSyntheticMeshBuilder::removeFaces(HelperTriangleMesh * mesh, size_t count)
    {
        std::vector<int>& indices = mesh->getIndices();

        for (size_t i = 0; i < count; i++)
        {
            size_t face = rand() % (indices.size() / 3);

            indices[face * 3 + 0] = indices[indices.size() - 3];
            indices[face * 3 + 1] = indices[indices.size() - 2];
            indices[face * 3 + 2] = indices[indices.size() - 1];
            indices.resize(indices.size() - 3);
        }
    }

    inline void HelperSyntheticMeshBuilder::addRandomness(std::vector<MathVector3f> & vertices, float randomness)
    {
        if (randomness == 0.0)
            return;

        for (size_t i = 0; i < vertices.size(); i++)
        {
            MathVector3f random((float)rand() / (float)(RAND_MAX), (float)rand() / (float)(RAND_MAX), (float)rand() / (float)(RAND_MAX));
            random.normalize();
            random *= randomness;
            vertices[i] += random;
        }
    }

    inline void HelperSyntheticMeshBuilder::rescaleToUnitBox(HelperTriangleMeshContainer * scene)
    {
        auto myMeshArray = scene->getMeshArray();

        std::vector<MathVector3f> myPoints;
        for (auto iter = myMeshArray.begin(); iter != myMeshArray.end(); iter++)
        {
            MathVector3f myMin, myMax;
            MathArithmetic<float>::getMinMax((*iter)->getVertices(), myMin, myMax);

            myPoints.push_back(myMin);
            myPoints.push_back(myMax);
        }

        MathVector3f myMin;
        MathVector3f myMax;

        MathArithmetic<float>::getMinMax(myPoints, myMin, myMax);

        MathVector3f myMiddle = (myMin + myMax) * (-0.5);
        MathVector3f myExtent = myMax - myMin;
        float m = myExtent.x;
        m = m > myExtent.y ? m : myExtent.y;
        m = m > myExtent.z ? m : myExtent.z;
        m = 1.0f / m;

        for (auto iter = myMeshArray.begin(); iter != myMeshArray.end(); iter++)
        {
            translate((*iter), myMiddle);
            scale((*iter), m);
        }
    }

    inline void HelperSyntheticMeshBuilder::scale(std::vector<MathVector3f> & vertices, float scale)
    {
        for (size_t i = 0; i < vertices.size(); i++)
        {
            vertices[i] *= scale;
        }
    }

    inline void HelperSyntheticMeshBuilder::scale(HelperTriangleMeshContainer * scene, float aScale)
    {
        auto myMeshArray = scene->getMeshArray();

        for (auto iter = myMeshArray.begin(); iter != myMeshArray.end(); iter++)
        {
            scale((*iter), aScale);
        }
    }

    inline void HelperSyntheticMeshBuilder::rotate(HelperTriangleMeshContainer * scene, float pitch, float yaw, float roll)
    {
        auto myMeshArray = scene->getMeshArray();

        for (auto iter = myMeshArray.begin(); iter != myMeshArray.end(); iter++)
        {
            rotate((*iter), pitch, yaw, roll);
        }
    }

    inline void HelperSyntheticMeshBuilder::translate(HelperTriangleMeshContainer * scene, const MathVector3f & translation)
    {
        auto myMeshArray = scene->getMeshArray();

        for (auto iter = myMeshArray.begin(); iter != myMeshArray.end(); iter++)
        {
            translate((*iter), translation);
        }
    }

    inline void HelperSyntheticMeshBuilder::rotate(std::vector<MathVector3f> & vertices, float pitch, float yaw, float roll)
    {
        MathMatrixf pitchMatrix;
        pitchMatrix.setRotateX(pitch);

        MathMatrixf yawMatrix;
        yawMatrix.setRotateY(yaw);

        MathMatrixf rollMatrix;
        rollMatrix.setRotateZ(roll);

        for (size_t i = 0; i < vertices.size(); i++)
        {
            //TODO replace this by a single matrix concatenation
            vertices[i] = pitchMatrix.multiply(vertices[i]);
            vertices[i] = yawMatrix.multiply(vertices[i]);
            vertices[i] = rollMatrix.multiply(vertices[i]);
        }
    }

    inline void HelperSyntheticMeshBuilder::translate(std::vector<MathVector3f> & vertices, const MathVector3f & translation)
    {
        for (size_t i = 0; i < vertices.size(); i++)
        {
            vertices[i] += translation;
        }
    }

    inline void HelperSyntheticMeshBuilder::addRandomness(HelperTriangleMesh * mesh, float randomness)
    {
        addRandomness(mesh->getVertices(), randomness);
    }

    inline void HelperSyntheticMeshBuilder::scale(HelperTriangleMesh * mesh, float s)
    {
        scale(mesh->getVertices(), s);
    }

    inline void HelperSyntheticMeshBuilder::rotate(HelperTriangleMesh * mesh, float pitch, float yaw, float roll)
    {
        rotate(mesh->getVertices(), pitch, yaw, roll);
    }

    inline void HelperSyntheticMeshBuilder::translate(HelperTriangleMesh * mesh, const MathVector3f & translation)
    {
        translate(mesh->getVertices(), translation);
    }
}
