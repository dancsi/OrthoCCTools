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

#ifndef GREEDYCOLORSORT_H
#define GREEDYCOLORSORT_H


#include "MaximumCliqueBase.h"


template<class Graph>
class GreedyColorSort {
public:
    typedef typename Graph::VertexSet VertexSet;
    typedef typename VertexSet::VertexId VertexId;
    typedef std::vector<VertexId> NumberedSet;
    typedef Graph GraphType;
    
protected:
    std::vector<VertexSet> colorSet;
    Graph* graph;
    Timers* timers;
    
public:
    GreedyColorSort() : graph(nullptr) {}
    
    void init(GraphType* g, Timers& t) {
        graph = g;
        timers = &t;
        size_t n = graph->getNumVertices();
        colorSet.resize(n);
        for (auto cc : colorSet)
            cc.reserve(n);
    }
    
    void assignVertexNumber(VertexSet& vs, NumberedSet& ns, size_t i, VertexId vert, VertexId num) {vs[i] = vert; ns[i] = num;}
    
    bool notEmpty(const NumberedSet& ns) const {return ns.size() > 0;}
    
    VertexId topNumber(const NumberedSet& ns) const {return ns.back();}
    
    VertexId topVertex(const NumberedSet& ns, const VertexSet& vs) const {return vs.back();}
    
    void popTop(NumberedSet& ns, VertexSet& vs) {ns.pop_back(); vs.pop();}
    
    // take clique c and candidate vertex set p as input
    // return numbered vertex set np as output
    void numberSort(const VertexSet& c, VertexSet& p, NumberedSet& np, unsigned int maxSize = 0) {
        size_t m = p.size();
        size_t numbers = 0;
            
        for (size_t i = 0; i < m; i++) {
            colorSet[i].clear();
            VertexId v = p[i];
            size_t k = 0;
            
            while (intersectionExists(v, colorSet[k])) 
                k++;
                
            colorSet[k].add(v);
            numbers = std::max(numbers, k+1);
        }
        
        np.resize(m);
        for (size_t k = 0, i = 0; k < numbers; k++) {
            for (size_t j = 0; j < colorSet[k].size(); j++) {
                VertexId v = colorSet[k][j];
                p[i] = v; 
                np[i++] = k+1;
            }
        }
    }

protected:    
    bool intersectionExists(VertexId p, const VertexSet& vertices) const {
        const VertexSet& neighbourhood = graph->adjacencyMatrix[p];
        size_t n = vertices.size();
        for (size_t i = 0; i < n; ++i) {
            if (neighbourhood[vertices[i]] == true)
                return true;
        }
        return false;
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(GreedyColorSort, "Greedy color sort");

#endif // GREEDYCOLORSORT_H
