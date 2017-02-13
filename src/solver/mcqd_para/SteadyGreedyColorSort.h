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

#ifndef STEADYGREEDYCOLORSORT_H
#define STEADYGREEDYCOLORSORT_H


#include "MaximumCliqueBase.h"


// Set of vertices is based on std::vector
// every set consists of two sets with different ordering
// used as VertexSetRepresentation in MaximumCliqueProblem
template<class T>
class SteadyVectorSet : public std::vector<std::pair<T,T>> {
public:
    typedef T VertexId;
    typedef std::vector<std::pair<T,T>> VecPairT;
    
    using VecPairT::resize;
    using VecPairT::reserve;
    void add(const T& value) {push_back(std::make_pair(value));}
    T pop() {T temp=this->back(); this->pop_back(); return temp;}
    void remove(T value) {
        if (this->back() == value) this->pop_back();
        else {
            auto temp = std::find(this->begin(), this->end(), value);
            if (temp != this->end())
                this->erase(temp);
        }
    }
    using VecPairT::size;
    //VertexId operator[] (size_t i) const {return VecPairT::operator[](i).second();}
    //VertexId& operator[] (size_t i) {return VecPairT::operator[](i).second();}
    using VecPairT::operator[];
    using VecPairT::clear;
    using VecPairT::pop_back;
    VertexId back() const {return VecPairT::back().second;}
    const std::pair<T,T>& backPair() const {return VecPairT::back();}
};
REGISTER_TEMPLATE1_CLASS_NAME(SteadyVectorSet, "Vector based set of pairs");


template<class Graph>
class SteadyGreedyColorSort {
public:
    typedef typename Graph::VertexSet VertexSet;
    typedef typename VertexSet::VertexId VertexId;
    typedef SteadyVectorSet<VertexId> NumberedSet;
    typedef Graph GraphType;
    
protected:
    std::vector<VertexSet> colorSet;
    Graph* graph;
    
public:
    SteadyGreedyColorSort() : graph(nullptr) {}
    
    void init(GraphType* g) {
        graph = g;
        size_t n = graph->getNumVertices();
        colorSet.resize(n);
        for (auto cc : colorSet)
            cc.reserve(n);
    }
    
    void assignVertexNumber(VertexSet& vs, NumberedSet& ns, size_t i, VertexId vert, VertexId num) {ns[i] = std::make_pair(vert, num);}
    
    bool notEmpty(const NumberedSet& ns) const {return ns.size() > 0;}
    
    VertexId topNumber(const NumberedSet& ns) const {return ns.back();}
    
    VertexId topVertex(const NumberedSet& ns, const VertexSet& vs) const {return ns.backPair().first;}
    
    void popTop(NumberedSet& ns, VertexSet& vs) {
        vs.remove(ns.backPair().first);
        ns.pop_back(); 
    }
    
    // take clique c and candidate vertex set p as input
    // return numbered vertex set np as output
    // difference compared to regular GreedyColorSort: p is here only used as input
    // p is input, np is output (p does not change)
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
                np[i].first = v; 
                //p[i] = v;
                np[i++].second = k+1;
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
REGISTER_TEMPLATE1_CLASS_NAME(SteadyGreedyColorSort, "Steady greedy color sort");

#endif // STEADYGREEDYCOLORSORT_H
