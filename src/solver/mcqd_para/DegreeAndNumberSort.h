/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef DEGREEANDNUMBERSORT_H_INCLUDED
#define DEGREEANDNUMBERSORT_H_INCLUDED


#include "MaximumCliqueBase.h"


template<class ColorSort>
struct DegreeAndNumberSort : ColorSort {
    typedef typename ColorSort::GraphType GraphType;
    typedef typename GraphType::VertexSet VertexSet;
    typedef typename ColorSort::NumberedSet NumberedSet;
    
    using ColorSort::numberSort;
    using ColorSort::init;
    using ColorSort::assignVertexNumber;
    
    void initialSort(VertexSet& c, VertexSet& vertices, NumberedSet& color) {
        typedef typename GraphType::VertexId VertexId;
        typedef std::pair<VertexId,VertexId> VerDeg;
        
        // sort by degree
        size_t n = vertices.size();
        size_t maxDegree = 0;
        std::vector<VerDeg> vertexAndDegree(n);
        for (size_t i = 0; i < n; ++i) {
            vertexAndDegree[i] = std::make_pair(i, this->graph->degrees[i]);
            maxDegree = std::max(maxDegree, (size_t)this->graph->degrees[i]);
        }
        std::stable_sort(vertexAndDegree.begin(), vertexAndDegree.end(), [](const VerDeg& a, const VerDeg& b){return (a.second > b.second);});
        
        // reorder vertices in graph
        std::vector<VertexId> vertexSet;
        vertexSet.reserve(n);
        for (const auto& vAd : vertexAndDegree)
            vertexSet.push_back(vAd.first);
        this->graph->orderVertices(vertexSet);
        
        // initial number sort
        color.resize(n);
        for (size_t i = 0; i < n; ++i) 
            assignVertexNumber(vertices, color, i, i, 0); // color is set to 0 as it does not matter here
        numberSort(c, vertices, color, 0);
    }
};
REGISTER_TEMPLATE_EXT_CLASS_NAME(DegreeAndNumberSort, "Degree/Color sort");


#endif // DEGREEANDNUMBERSORT_H_INCLUDED
