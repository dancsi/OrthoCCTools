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

#ifndef BB_GREEDYCOLORSORT_H_INCLUDED
#define BB_GREEDYCOLORSORT_H_INCLUDED


#include "BitSet.h"
#include "MaximumCliqueBase.h"
#include "SteadyGreedyColorSort.h" // SteadyVectorSet is reused here


// Set of vertices is based on std::bitset
// every set consists of two sets with different ordering
// used as VertexSetRepresentation in MaximumCliqueProblem
class BitstringSet : protected BitSet {
    mutable size_t countCache;
    
public:
    typedef unsigned int VertexId;
    
    ~BitstringSet() {}
    BitstringSet() : countCache(0) {}
    using BitSet::resize;
    void reserve(size_t s) {BitSet::resize(s);};
    void add(VertexId value) {setValue(value, true); ++countCache;}
    void remove(VertexId value) {setValue(value, false); --countCache;}
    size_t recount() const {countCache = count(); return countCache;}
    size_t size() const {return countCache;}
    size_t alsize() const {return BitSet::size();}
    bool operator[](size_t index) const {return getValue(index);}
    BoolProxy operator[](size_t index)  {return BoolProxy(this, index);}
    //using BitSet::operator[];
    void clear() {BitSet::clear(); countCache = 0;};
    using BitSet::nextSetBit;
    using BitSet::operator~;
    
    void operator&=(const BitSet& other) {BitSet::operator&= (other);}
    void operator&=(const BitstringSet& other) {BitSet::operator&= (other);}
    void operator^=(const BitstringSet& other) {BitSet::operator^= (other);}
    
    friend void intersectWithAdjecency (const BitstringSet& v, const BitstringSet& adj, BitstringSet& result) {
        result = v;
        result &= adj;
        result.recount();
    }
    
    void remove(const BitstringSet& values) {
        // this function only works correctly if all the specified values are set in the target bitstring
        *this ^= values;
        countCache -= values.countCache;
    }
    
    bool isIntersectionOf(const BitstringSet& bigSet) {
        return  false;
    }
};
REGISTER_CLASS_NAME(BitstringSet, "bitstring based set");

template<> struct isTypeASet<BitstringSet> { enum {value = true}; };

template<class Out>
Out& operator<< (Out& out, const BitstringSet& b) {
    auto bCopy = b;
    if (b.size() == 0) {
        std::cout << "[/]";
    } else {
        std::cout << "[";
        bool comma = false;
        while (bCopy.size() > 0) {
            if (comma)
                std::cout << "," << std::flush;
            comma = true;
            auto v = bCopy.nextSetBit();
            std::cout << v;
            bCopy.remove(v);
        }
        std::cout << "]";
    }
    return out;
}


template<class Graph>
class BBGreedyColorSort {
public:
    typedef BitstringSet VertexSet;
    typedef typename VertexSet::VertexId VertexId;
    typedef SteadyVectorSet<VertexId> NumberedSet;
    typedef Graph GraphType;
    
protected:
    std::vector<VertexSet> colorSet;
    Graph* graph;
    VertexSet Ubb, Qbb;
    
public:
    BBGreedyColorSort() : graph(nullptr) {}
    
    void init(GraphType* g) {
        graph = g;
        
        size_t n = graph->getNumVertices();
        colorSet.resize(n);
        for (auto& cc : colorSet) {
            cc.reserve(n);
        }
    }
    
    void assignVertexNumber(VertexSet& vs, NumberedSet& ns, size_t i, VertexId vert, VertexId num) {ns[i] = std::make_pair(vert, num);}
    
    bool notEmpty(const NumberedSet& ns) const {return ns.size() > 0;}
    
    VertexId topNumber(const NumberedSet& ns) const {return ns.back();}
    
    VertexId topVertex(const NumberedSet& ns, const VertexSet& vs) const {return ns.backPair().first;}
    
    void popTop(NumberedSet& ns, VertexSet& vs) {auto v = ns.backPair().first; ns.pop_back(); vs.remove(v);}
    
    // take clique c and candidate vertex set p as input
    // return numbered vertex set np as output
    // difference compared to regular GreedyColorSort: p is here only used as input
    void numberSort(const VertexSet& c, const VertexSet& p, NumberedSet& np, unsigned int maxSize = 0) {
        Ubb = p;
        Qbb = Ubb;
        int k = 0;
        size_t i = 0;
        int kMin = (int)maxSize - (int)c.size();
        np.resize(p.size());
        
        while (Ubb.size() > 0) {
            while (Qbb.size() > 0) {
                auto v = Qbb.nextSetBit();
                Qbb &= graph->invAdjacencyMatrix[v];
                Qbb.recount();
                Ubb.remove(v);
                if (k >= kMin) {
                    np[i].first = v;
                    np[i++].second = k+1;
                }
            }
            Qbb = Ubb;
            ++k;
        }
        np.resize(i);
    }

protected:    
    bool intersectionExists(VertexId p, const VertexSet& vertices) const {
        const VertexSet& neighbourhood = graph->adjacencyMatrix[p];
        auto intersection = neighbourhood;
        intersection &= vertices;
        return intersection.size() > 0;
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(BBGreedyColorSort, "greedy color sort on bitstrings");

#endif // BB_GREEDYCOLORSORT_H_INCLUDED
