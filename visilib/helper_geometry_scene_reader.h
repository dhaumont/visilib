/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2021 by Denis Haumont

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
#include <fstream>
#include <iostream>
#include <string>

#include "helper_triangle_mesh_container.h"
#include "helper_triangle_mesh.h"

namespace visilib
{
    class HelperTriangleMeshContainer;

    /** @brief An helper class allowing to read triangle mesh stored in the ".obj" file format 
 
    The obj reader code is inspired from Meshlab implementation http://meshlab.sourceforge.net/
    */

    class HelperGeometrySceneReader
    {
    public:

        HelperGeometrySceneReader(HelperTriangleMeshContainer* scene);

        ~HelperGeometrySceneReader();

        bool readFileObj(const std::string& fileName);

        static void tokenizeNextLine(std::istream& stream, std::vector< std::string >& tokens);


    private:

        HelperTriangleMeshContainer* mScene;
    };

    inline HelperGeometrySceneReader::HelperGeometrySceneReader(HelperTriangleMeshContainer* aScene)
    {
        mScene = aScene;
    }

    inline HelperGeometrySceneReader::~HelperGeometrySceneReader()
    {
    }

    //--------------------------------------------------------------------------------------
    // obj reader code from Meshlab
    //http://meshlab.sourceforge.net/
    //--------------------------------------------------------------------------------------

    inline void HelperGeometrySceneReader::tokenizeNextLine(std::istream& stream, std::vector< std::string >& tokens)
    {
        if (stream.eof())
            return;

        std::string line;
        do
            std::getline(stream, line);
        while ((line.length() == 0 || line[0] == '#') && !stream.eof());  // skip comments and empty lines

        if ((line[0] == '#') || (line.length() == 0))  // can be true only on last line of file
            return;

        int from = 0;
        int to = 0;
        size_t length = line.size();
        tokens.clear();
        do
        {
            while (from != length && (line[from] == ' ' || line[from] == '\t' || line[from] == '\r'))
                from++;
            if (from != length)
            {
                to = from + 1;
                while (to != length && line[to] != ' ' && line[to] != '\t' && line[to] != '\r')
                    to++;
                tokens.push_back(line.substr(from, to - from).c_str());
                from = to;
            }
        } while (from < length);
    } // end TokenizeNextLine

    inline bool HelperGeometrySceneReader::readFileObj(const std::string & fileName)
    {
        std::ifstream stream(fileName.c_str());

        if (!stream.good())
            return false;

        std::vector< std::string > tokens;
        std::string	header;

        std::vector<float> vertices;
        std::vector<int> indices;

        size_t faceCount = 0;
        size_t firstVertexIndex = 0;
        size_t indexCount = 0;

        while (!stream.eof())
        {
            tokens.clear();
            tokenizeNextLine(stream, tokens);

            int numTokens = static_cast<int>(tokens.size());
            if (numTokens == 0) continue;

            header.clear();
            header = tokens[0];

            if (header.compare("g") == 0)	// vertex
            {
                /*std::string surfaceName = "";
                for (int i = 1; i < numTokens; i++)
                {
                surfaceName += tokens[i];
                surfaceName += " ";
                }

                if (faceCount > 0 && vertices.size() > 0)
                {
                GeometryMesh* mesh = new GeometryMesh(vertices, indices);

                mScene->add(mesh);
                firstVertexIndex += (vertices.size()/3);
                }

                faceCount = 0;
                indexCount = 0;

                vertices.clear();
                indices.clear();			*/
            }
            if (header.compare("v") == 0)	// vertex
            {
                if (numTokens < 4) return false;

                float x = (float)atof(tokens[1].c_str());
                float y = (float)atof(tokens[2].c_str());
                float z = (float)atof(tokens[3].c_str());
                vertices.push_back(x);
                vertices.push_back(y);
                vertices.push_back(z);
            }
            else if (header.compare("vt") == 0)	// vertex texture coords
            {
                if (numTokens < 3) return false;

                //float u = static_cast<float>(atof(tokens[1].c_str()));
                //float t = static_cast<float>(atof(tokens[2].c_str()));
            }
            else if (header.compare("vn") == 0)  // vertex normal
            {
                if (numTokens != 4) return false;

                //float nx = (float)atof(tokens[1].c_str());
                //float ny = (float)atof(tokens[2].c_str());
                //float nz = (float)atof(tokens[3].c_str());
            }
            else if ((header.compare("f") == 0))  // face
            {
                size_t vertexCount = tokens.size() - 1;

                if (vertexCount >= 3)
                {
                    faceCount++;

                    for (int i = 0; i < vertexCount; ++i) // remember index starts from 1 instead of 0
                    {
                        size_t vertexIndex = atoi(tokens[i + 1].c_str());
                        indices.push_back(static_cast<int>(vertexIndex - 1 - firstVertexIndex));
                    }
                    indexCount += vertexCount;
                }
            }
        }

        if (faceCount > 0 && vertices.size() > 0)
        {
            if (indices[0] < 0)
            {
                size_t firstVertex = vertices.size() / 3 + 1;

                for (int i = 0; i < indices.size(); i++)
                {
                    indices[i] = static_cast<int>(firstVertex + indices[i]);
                }
            }
            HelperTriangleMesh* mesh = new HelperTriangleMesh(vertices, indices);

            mScene->add(mesh);
        }

        return true;
    }
};
