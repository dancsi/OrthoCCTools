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

#ifndef HEADER_38B6D205806EECF7
#define HEADER_38B6D205806EECF7


#include "BB_GreedyColorSort.h"


template<class Graph>
class BBColorRSort : public BBGreedyColorSort<Graph> {
public:
    typedef BBGreedyColorSort<Graph> ParentT;
    typedef typename ParentT::VertexSet VertexSet;
    typedef typename ParentT::NumberedSet NumberedSet;
    typedef typename ParentT::VertexId VertexId;
    
protected:
    // tell compiler about parent class member variables
    using ParentT::Ubb;
    using ParentT::Qbb;
    VertexSet intersection; // used in resort function
    
public:
    using ParentT::colorSet;
    using ParentT::graph;
    
    void numberSort(const VertexSet& c, const VertexSet& p, NumberedSet& np, unsigned int maxSize = 0) {
        Ubb = p;
        Qbb = Ubb;
        int k = 0;
        size_t i = 0;
        int kMin = (int)maxSize - (int)c.size();
        np.resize(p.size());
        
        while (Ubb.size() > 0) {
            bool first = false;
            colorSet[k].clear();
            Qbb = Ubb;
            while (Qbb.size() > 0) {
                auto v = Qbb.nextSetBit();
                if ((Qbb.size() == 1) && (k >= kMin) && first) {
                    //std::cout << "Recolor " << v << " of Qbb " << Qbb.size() << " " << k << " " << kMin << "\n";
                    if (recolor(v, Qbb, k, kMin)) {
                        Ubb.remove(v);
                        // because of the condition (Qbb.cachedSize() == 1), Qbb only contains v at this point
                        // so jump right out of the loop to save several cycles
                        break;
                        //Qbb.remove(v);
                    } else {
                        first = false;
                    }
                }
                colorSet[k].add(v);
                Qbb &= graph->invAdjacencyMatrix[v];
                Qbb.recount(); 
            }
            Ubb.remove(colorSet[k]);
            ++k;
        }
        
        for (int k1 = std::max(0,kMin); k1 < k; ++k1) {
            for (int j = colorSet[k1].size(); j > 0; --j) {
                auto v = colorSet[k1].nextSetBit();
                colorSet[k1].remove(v);
                np[i].first = v;
                np[i++].second = k1+1;
            }
        }
        np.resize(i);
    }
    
    bool recolor(VertexId v, VertexSet& p, int kv, int kMin) {
        for (int k1 = 0; k1 < kMin - 1; ++k1) {
            intersection = colorSet[k1];
            intersection &= graph->adjacencyMatrix[v];
            intersection.recount();
            if (intersection.size() == 1) {
                auto w = intersection.nextSetBit();
                for (int k2 = k1+1; k2 < kMin; ++k2) {
                    intersection = colorSet[k2];
                    intersection &= graph->adjacencyMatrix[w];
                    intersection.recount();
                    if (intersection.size() == 0) {
                        colorSet[k1].add(v);
                        colorSet[k1].remove(w);
                        colorSet[k2].add(w);
                        return true;
                    }
                }
            } else if (intersection.size() == 0) {
                colorSet[k1].add(v);
                return true;
            }
        }
        return false;
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(BBColorRSort, "greedy color sort with resort on bitstrings");


#endif // header guard 
