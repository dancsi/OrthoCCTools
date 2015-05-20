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

#ifndef MCRBB_H_INCLUDED
#define MCRBB_H_INCLUDED

/**
    This file is a copy of Mcr.h but fixed for BB. TODO: include support for BBs to Mcr
**/

#include "MaximumCliqueBase.h"


template<class ColorSort>
struct BBMcrSort : public ColorSort {
    using ColorSort::numberSort;
    typedef typename ColorSort::GraphType GraphType;
    typedef typename GraphType::VertexSet VertexSet;
    typedef typename VertexSet::VertexId VertexId;
    typedef typename ColorSort::NumberedSet NumberedSet;
    
    struct Vertex {
        int degree;
        int exDegree;
        int index;
        template<class Ostream>
        friend Ostream& operator<< (Ostream& out, const Vertex& v) {
            out << "<" << v.index << ":" << v.degree << ";" << v.exDegree << ">";
            return out;
        }
    };
    
    using ColorSort::init;
    using ColorSort::assignVertexNumber;
    
    BBMcrSort() {}
        
    void initialSort(VertexSet& c, VertexSet& vertices, NumberedSet& color) {
        TRACE("MCR style initial sort", TRACE_MASK_CLIQUE | TRACE_MASK_INITIAL, 1); 
        size_t n = vertices.size();
        if (n == 0) {
            TRACE("initial sort foud n == 0 and exited", TRACE_MASK_CLIQUE | TRACE_MASK_INITIAL, 1); 
            return;
        }
        auto maxDegree = this->graph->degrees[0];
        std::vector<Vertex> r(n);
        std::vector<VertexId> order(n); // reordering vector
        
        for (size_t i = 0; i < n; ++i) {
            r[i].index = i;
            r[i].degree = this->graph->degrees[r[i].index];
            r[i].exDegree = 0;
            maxDegree = std::max(maxDegree, r[i].degree);
        }
        TRACEVAR(r, TRACE_MASK_INITIAL, 2);
        
        // not sure if the following calculation is correct for ex-deg (not clearly specified in Tomita 2006)
        //  it is possible this should be done on every step of the following while loop, taking only
        //  the neighbourhood of the observed vertex into an account (but probably not)
        for (size_t i = 0; i < n; ++i) 
            for (size_t j = 0; j < n; ++j) 
                if (this->graph->adjacencyMatrix[i][j] == true)
                    r[i].exDegree += r[j].degree;
                    
        // sort by degree (descending), stable mode (respect relative order of the vertices with the same degree)
        std::stable_sort(r.begin(), r.end(), [](const Vertex& a, const Vertex& b) {return (a.degree > b.degree);});
        TRACEVAR(r, TRACE_MASK_INITIAL, 2); 
                
        // index in vertices
        size_t vi = n-1;
        
        size_t rMinIndex = r.size()-1;
        while (rMinIndex > 0) {
            // locate vertices with min degree
            // set of vertices "Rmin" is implemented as a subarray from index rMinIndex to the end of the set r
            int minDeg = r.back().degree;
            rMinIndex = r.size()-1;
            while ((rMinIndex > 0) && (r[rMinIndex-1].degree == minDeg)) 
                --rMinIndex;
            
            if (rMinIndex == 0)
                break;
                
            // if "Rmin" contains more than 1 element
            if (rMinIndex < r.size()-1) {
                // sort by ex-deg (descending - max is first)
                std::stable_sort(r.begin()+rMinIndex, r.end(), [](const Vertex& a, const Vertex& b){return (a.exDegree > b.exDegree);});
            }
            
            // vertex with min ex-deg in rMin goes into the ordered set of vertices (filled from the back towards the front)
            Vertex& p = r.back();
            order[vi] = p.index;
            --vi;
            
            // decrease the degree of remaining vertices that are adjacent to p
            auto& adjacent = this->graph->adjacencyMatrix[p.index];
            r.pop_back();
            size_t rs = r.size();
            for (size_t i = rs; i > 0; --i) {
                if (adjacent[r[i-1].index] == true) {
                    auto rid = --r[i-1].degree;
                    // sort the modified vertex immediately
                    for (size_t j = i; (j < rs) && (rid  < r[j].degree); ++j) 
                        swap(r[j-1], r[j]);
                }
            }
        }
        TRACE("degree of leftover vortices MCR initial sort:", TRACE_MASK_INITIAL, 2); 
        TRACEVAR(r[0].degree, TRACE_MASK_INITIAL, 2);
        TRACE("after calculation of ex-degree:", TRACE_MASK_INITIAL, 2); 
        TRACEVAR(r, TRACE_MASK_INITIAL, 2); 
        TRACEVAR(c, TRACE_MASK_INITIAL, 2); 
//      std::cout << "DEBUG: iterated to regular subgraph of degree " << r.front().degree << " and size " << vi << "\n";
        
        // all the vertices in r have the same degree (regular subgraph) → perform ordinary number sort (color sort) 
        c.clear();
        c.reserve(n);
        VertexSet dummySet;
        dummySet.reserve(n);
        
        for (size_t i = 0; i < r.size(); ++i) {
            c.add(r[i].index);
            dummySet.add(r[i].index);
        }
        
        numberSort(dummySet, c, color, 0);
        
        // calculate maximum number (color), then fill in the rest of the order vector (clearing colorC and c in the process)
        size_t m = color.back();
        size_t mmax = r.size() + maxDegree - m;
        std::vector<VertexId> storedColor(r.size());
        for (size_t i = r.size(); i > 0; --i) {
            order[i-1] = this->topVertex(color, c); 
            storedColor[i-1] = this->topNumber(color); 
            this->popTop(color, c);
        }

        // if the degree of vertices in r (which all have the same degree) equals r.size() - 1, then r is a clique
        if (r.size() == r[0].degree+1) {
        	// fill the clique with the new vertex numbers (these are actually numbers [0..r.size()-1]), the ones that will be set by the orderVertices function also for the rest of the graph
        	//c.resize(r.size());
        	c.reserve(n);
        	for (size_t i = 0; i < r.size(); ++i) {
        		c.add(i);
        	}
        	TRACE("initial clique found:", TRACE_MASK_INITIAL, 2); 
        	TRACEVAR(c, TRACE_MASK_INITIAL, 2); 
        }
        
        // reorder vertices in graph
        this->graph->orderVertices(order);
        
        // number
        color.resize(n); 
        vertices.clear();
        // first few colors remain as they are, the rest are filled in with up to mmax
        for (size_t i = 0; i < n; ++i) {
            if (i < r.size()) {
                assignVertexNumber(vertices, color, i, i, storedColor[i]);
            } else if (i < mmax) {
                ++m;
                assignVertexNumber(vertices, color, i, i, m);
            } else {
                assignVertexNumber(vertices, color, i, i, maxDegree + 1);
            }
            vertices.add(i);
            //std::cout << i << ": " << this->graph->mapping[color[i].first] << "," << color[i].second << "  ";
        }
        
        TRACE("initial coloring of vertices:", TRACE_MASK_INITIAL, 2); 
        TRACEVAR(color, TRACE_MASK_INITIAL, 2); 
        TRACEVAR(c, TRACE_MASK_INITIAL, 2); 
    }
};
REGISTER_TEMPLATE_EXT_CLASS_NAME(BBMcrSort, "MCR sort for Bitstrings");


#endif // MCRBB_H_INCLUDED
