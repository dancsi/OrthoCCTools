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

#ifndef MAXIMUMCLIQUEBASE_H
#define MAXIMUMCLIQUEBASE_H


#include "KillTimer.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <bitset>

#include "Timer.h"

template<class CharPtr>
CharPtr skipPath(CharPtr text) {
    CharPtr lastSlash = text, t = text;
    while (*t != 0) {
        if (*t == '\\' || *t == '/') lastSlash = t+1;
        ++t; 
    }
    return lastSlash;
}

#ifdef NDEBUG1
    #define TRACE(text, mask, level) {}
    #define TRACEVAR(var, mask, level) {}
#else
    #define IF_TRACABLE(mask, level) if ((level <= TRACE_LEVEL) && ((TRACE_MASK & mask) > 0))
    #define TRACE(text, mask, level) {IF_TRACABLE(mask, level) std::cout  << "(" << level << ") " << skipPath(__FILE__) << " line " << __LINE__ << ": " << text << "\n";}
    #define TRACEVAR(var, mask, level) {IF_TRACABLE(mask, level) std::cout  << "(" << level << ") " << skipPath(__FILE__) << " line " << __LINE__ << ": " << #var << "=" << var << "\n";}
#endif

    
#define TRACE_MASK_THREAD 0x10
#define TRACE_MASK_CLIQUE 0x20
#define TRACE_MASK_INITIAL 0x40



template<class Ostream, class T>
Ostream& operator<< (Ostream& out, const std::vector<T>& vec) {
    if (vec.size() > 0) {
        out << "[" << vec[0];
        for (size_t i = 1; i < vec.size(); ++i)
            out << "," << vec[i]; 
        out << "]";
    } else 
        out << "[/]";
    return out;
}


template<class Ostream, class T, class U>
Ostream& operator<< (Ostream& out, const std::vector<std::pair<T, U>>& vec) {
    if (vec.size() > 0) {
        out << "[<" << vec[0].first << "," << vec[0].second << ">";
        for (size_t i = 1; i < vec.size(); ++i)
            out << "; <" << vec[i].first << "," << vec[i].second << ">"; 
        out << "]";
    } else 
        out << "[/]";
    return out;
}

template<class T>
struct ClassName {
    static const char* getValue() {
        return "";
    }
};

#define REGISTER_CLASS_NAME(CLASS, NAME) template<> struct ClassName<CLASS> {static const char* getValue() {return NAME;} };
#define REGISTER_TEMPLATE1_CLASS_NAME(CLASS, NAME) template<class T> struct ClassName<CLASS<T>> {static const char* getValue() {return NAME;} };
#define REGISTER_TEMPLATE2_CLASS_NAME(CLASS, NAME) template<class T1, classT2> struct ClassName<CLASS<T1,T2>> {static const char* getValue() {return NAME;} };
#define REGISTER_TEMPLATE_EXT_CLASS_NAME(CLASS, NAME) template<class T> struct ClassName<CLASS<T>> {static const char* getValue() { \
    static std::unique_ptr<std::string> rs(new std::string); \
    if (rs->empty()) { std::stringstream s; s << NAME << "<" << ClassName<T>::getValue() << ">"; *rs = s.str(); } return rs->c_str(); } };


template<
    class VertexRepresentation,
    class SetRepresentation, 
    class GraphRepresentation, 
    class Sorter
>
class MaximumCliqueProblem : Sorter {
public:
    typedef VertexRepresentation VertexId;
    typedef SetRepresentation VertexSet;
    
protected:
    using Sorter::notEmpty;
    using Sorter::topNumber;
    using Sorter::topVertex;
    using Sorter::popTop;
    using Sorter::numberSort;
    using Sorter::initialSort;
    typedef typename Sorter::NumberedSet NumberedSet;
        
    std::string algorithmName;
    GraphRepresentation* graph;
    VertexId n;                         // number of vertices
    unsigned int maxSize;               // size of max clique
    VertexSet maxClique;
    PrecisionTimer timer;
    unsigned long long int steps;
    KillTimer1 killTimer;
    
public:
    VertexSet knownC;

    MaximumCliqueProblem(GraphRepresentation& graph) : graph(&graph), n(graph.getNumVertices()), maxSize(0), steps(0) {}
    
    // get the result of the search - maximal clique of the provided graph
    const VertexSet& getClique() const {return maxClique;}
    
    // output statistics of the last search (mostly timer readings); colored=true produces colored text on terminals supporting ANSI escape sequences
    void outputStatistics(bool colored = true) {
        std::ostringstream algorithmName;
        auto basefmt = std::cout.flags();
        auto baseFill = std::cout.fill();
        algorithmName << "MC(" << ClassName<VertexId>::getValue() << ","
            << ClassName<VertexSet>::getValue() << "," 
            << ClassName<GraphRepresentation>::getValue() << "," << ClassName<Sorter>::getValue() << ") ";
        std::cout << "-- " << std::setw(80-3) << std::setfill('-') << std::left << algorithmName.str();
        std::cout.flags(basefmt);
        std::cout << std::setfill(baseFill) << "\n";
        std::cout << "search took " << timer.totalSeconds() << "s; " << steps << " steps\n";
        VertexSet c = getClique();
        graph->remap(c);
        if (wasSearchInterrupted()) std::cout << "Warning, search has been interrupted, the results might not be correct\n";
        std::cout << "Clique (" << getClique().size() << "): " << c;
    }
    
    // run the search for max clique
    void search() {
        killTimer.start(10000);
        ScopeTimer t(timer);
        VertexSet c; // clique
        VertexSet p; // working set of vertices
        NumberedSet numbers;       // ordered set of colors
        
        c.reserve(n);
        p.reserve(n);
        for (VertexId i = 0; i < n; i++) 
            p.add(i);
            
        Sorter::init(graph);
        initialSort(c, p, numbers);
        
        TRACEVAR(c, TRACE_MASK_INITIAL, 2); 

        // some initial sorts (e.g. MCR) also find a clique
        maxSize = c.size();
        if (maxSize > 0)
            saveSolution(c);
                            
        if (numbers.size() == 0)
            numberSort(c, p, numbers, maxSize);
        TRACEVAR(numbers, TRACE_MASK_INITIAL, 2); 
        c.clear();
        
        expand(c, p, numbers);
        
        killTimer.cancel();
    }
    
    bool wasSearchInterrupted() const {return killTimer.timedOut;}
    
protected:   
    // main recursive function 
    //  c ... clique candidate, 
    //  p ... set of candidate vertices that can be added to c,
    //  numbered ... an ordered and numbered set of vertices
    void expand(VertexSet& c, VertexSet& p, NumberedSet& np) {
        ++steps;
        while (notEmpty(np)) {
            if (c.size() + topNumber(np) <= maxSize || killTimer.timedOut) {return;}
            auto v = topVertex(np, p);
            //std::cout << "v" << v << " ";
            c.add(v);
            VertexSet p1;
            graph->intersectWithNeighbours(v, p, p1);
            
            if (p1.size() == 0) {
                if (c.size() > maxSize) 
                    saveSolution(c);
            } else {
                NumberedSet np1;
                numberSort(c, p1, np1, maxSize);
                expand(c, p1, np1);
            }
            
            c.remove(v);
            popTop(np, p);
            //std::cout << "p" << np.size() << " ";
        }
    }
    
    // when a clique, larger than its predecessor is found, call this function to store it
    void saveSolution(const VertexSet& c) {
        maxSize = c.size();
        maxClique = c;
        //std::cout << "new clique " << maxSize << "\n";
    }
};

template<class T>
struct isTypeASet {
    enum {value = false};
};

template<class VectorSetRepresentation>
class Graph {
public:
    std::vector<VectorSetRepresentation> adjacencyMatrix, invAdjacencyMatrix;
    std::vector<int> degrees;
    std::vector<int> mapping;
    
    typedef typename VectorSetRepresentation::VertexId VertexId;
    typedef VectorSetRepresentation VertexSet;
    
    size_t getNumVertices() const {return adjacencyMatrix.size();}
    void init(const std::vector<std::vector<char> >& adjacency, const std::vector<int>& d) {
        size_t n = adjacency.size();
        adjacencyMatrix.resize(n); 
        invAdjacencyMatrix.resize(n); 
        for (size_t i = 0; i < n; ++i) {
            adjacencyMatrix[i].resize(n, false);
            invAdjacencyMatrix[i].resize(n, false);
            invAdjacencyMatrix[i][i] = false;
            for (size_t j = i+1; j < adjacency[i].size(); ++j) {
                adjacencyMatrix[i][j] = adjacency[i][j];
                invAdjacencyMatrix[i][j] = (adjacency[i][j] == false);
            }
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                adjacencyMatrix[i][j] = adjacency[j][i];
                invAdjacencyMatrix[i][j] = (adjacency[j][i] == false);
            }
        }
        degrees = d;
        mapping.clear();
    }
    
    // perform intersection, "result" keeps the ordering of the set "vertices"
    void intersectWithNeighbours(VertexId p, const VertexSet& vertices, VertexSet& result) const {
        // global function intersectWithAdjecency(VectorSetRepresentation, VectorSetRepresentation, VectorSetRepresentation) must be specified
        intersectWithAdjecency(vertices, adjacencyMatrix[p], result);
    }
    
    bool intersectionExists(VertexId p, const VertexSet& vertices) const {
        size_t n = vertices.size();
        for (size_t i = 0; i < n; ++i) {
            if (adjacencyMatrix[p][vertices[i]])
                return true;
        }
        return false;
    }
    
    // change the order of vertices in the adjacency matrix (renumber them)
    template<class Vec>
    void orderVertices(const Vec& order) {
        // check order vector
        size_t n = getNumVertices();
        if (order.size() != n) 
            throw "Invalid size vector in orderVertices";
            
        // create a vertex mapping table that will be used to renumber vertices back to original
        if (mapping.size() == 0) {
            // create default mapping that maps i → i
            mapping.resize(n);
            for (size_t i = 0; i < n; ++i)
                mapping[i] = i;
        }
        decltype(mapping) mapping2(n);
        
        // remap to temporary adjacencyMatrix
        std::vector<VectorSetRepresentation> adjacencyMatrix2;
        adjacencyMatrix2.resize(n);
        for (size_t i = 0; i < n; ++i) {
            adjacencyMatrix2[i].resize(n);
            invAdjacencyMatrix[i].clear();
            invAdjacencyMatrix[i].resize(n);
            mapping2[i] = mapping[order[i]];
            auto& adjRowI = adjacencyMatrix[order[i]];
            for (size_t j = 0; j < n; ++j) {
                adjacencyMatrix2[i][j] = adjRowI[order[j]] == true;
                invAdjacencyMatrix[i][j] = (i != j) & !adjacencyMatrix2[i][j];
            }
            // the line above includes the condition (i != j) because:
            // adjacency inverse is used to filter out vertices (operator &) and it is useful if 
            // given a vertex, it filters out its neighbours as well as the vertex itself
            // therefore reset the edge linking vertices to themselves
        }
        std::swap(adjacencyMatrix2, adjacencyMatrix);
        std::swap(mapping2, mapping);
        auto oldDeg = degrees;
        for (size_t i = 0; i < degrees.size(); ++i) degrees[i] = oldDeg[order[i]];
    }
    
    // in-place remapping function
    void remap(VectorSetRepresentation& v) {
        if (mapping.size() == 0 || v.size() == 0) return;
        
        if (isTypeASet<VectorSetRepresentation>::value) {
            VectorSetRepresentation rv;
            rv.reserve(mapping.size());
            
            if (mapping.size() < v.size())
                throw "Mapping failed, mapping is not known for all vertices";
            
            for (size_t i = 0; rv.size() < v.size(); ++i) {
                if (v[i]) rv.add(mapping[i]);
            }
            std::swap(rv, v);
        } else
            for (size_t i = 0; i < v.size(); ++i) v[i] = mapping[v[i]];
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(Graph, "Graph");

template<class G, class S>
bool isClique(const G& graph, const S& v) {
    size_t n = v.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (graph.adjacencyMatrix[v[i]][v[j]] == false) {
                //std::cout << v[i] << " and " << v[j] << " are not connected\n" << graph.adjacencyMatrix[v[i]][v[j]] << "\n";
                return false;
            }
        }
    }
    return true;
}

// Set of vertices is based on std::vector
// used as VertexSetRepresentation in MaximumCliqueProblem
template<class T>
class VectorSet : std::vector<T> {
public:
    typedef T VertexId;
    
    using std::vector<T>::resize;
    using std::vector<T>::reserve;
    void add(const T& value) {push_back(value);}
    T pop() {T temp=this->back(); this->pop_back(); return temp;}
    void remove(const T& value) {
        if (this->back() == value) this->pop_back();
        else {
            auto temp = std::find(this->begin(), this->end(), value);
            if (temp != this->end())
                this->erase(temp);
        }
    }
    using std::vector<T>::size;
    using std::vector<T>::operator[];
    using std::vector<T>::clear;
    using std::vector<T>::back;
    template<class AdjSet>
    friend void intersectWithAdjecency (const VectorSet& v, const AdjSet& adj, VectorSet& result) {
        auto n = v.size();
        result.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            if (adj[v[i]])
                result.add(v[i]);
        }
    }

    // only required for debugging
    bool isIntersectionOf(const VectorSet& bigSet) {
        size_t n = bigSet.size();
        if (n < 1) return false;
        --n;
        for (size_t i = 0; i < size(); ++i) {
            for (size_t j = 0; bigSet[j] != (*this)[i]; ++j) {
                if (j == n) return false;
            }
        }
        return true;
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(VectorSet, "Vector based set");

template<class Ostream, class T>
Ostream& operator<< (Ostream& out, const VectorSet<T>& vec) {
    if (vec.size() > 0) {
        out << "[" << vec[0];
        for (size_t i = 1; i < vec.size(); ++i)
            out << "," << vec[i]; 
        out << "]";
    } else 
        out << "[/]";
    return out;
}

REGISTER_CLASS_NAME(int, "int");
REGISTER_CLASS_NAME(short int, "short int");
REGISTER_CLASS_NAME(char, "char");

#endif // MAXIMUMCLIQUEBASE_H
