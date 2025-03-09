//
// Created by Yuzhen Feng on 19/2/2025.
//

#ifndef DEMAND_H
#define DEMAND_H

#include <vector>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "graph.h"

template<typename cost_type>
struct Path {
    std::list<typename Graph<cost_type>::vertex_type> vertex_list;
    std::list<typename Graph<cost_type>::edge_type> edge_list;
    boost::numeric::ublas::compressed_matrix<bool> edge_matrix;
    double flow;
    double cost;

    explicit Path(const Graph<cost_type>& graph) : flow(0.0), cost(0) {
        edge_matrix.resize(graph.num_vertices, graph.num_vertices, false);
    }

    void initialize(const Graph<cost_type>& graph,
                    const std::vector<typename Graph<cost_type>::vertex_type>& predecessors,
                    typename Graph<cost_type>::vertex_type origin,
                    typename Graph<cost_type>::vertex_type destination) {
        for (auto v = destination; v != origin; v = predecessors[v]) {
            vertex_list.push_front(v);
            edge_list.push_front(boost::edge(predecessors[v], v, graph.g).first);
            edge_matrix(static_cast<int>(predecessors[v]), static_cast<int>(v)) = true;
        }
        vertex_list.push_front(origin);
    }

    bool operator==(const Path& other) const {
        bool is_same = true;
        auto it1 = vertex_list.begin();
        auto it2 = other.vertex_list.begin();
        while (it1 != vertex_list.end() && it2 != other.vertex_list.end()) {
            if (*it1 != *it2) {
                is_same = false;
                break;
            }
            ++it1;
            ++it2;
        }
        return is_same;
    }

    void update_cost(const Graph<cost_type>& graph) {
        cost = 0.0;
        for (auto& edge : edge_list) {
            cost += graph.g[edge].cost;
        }
    }

    void print_path() {
        for (auto& v : vertex_list) {
            std::cout << v << " > ";
        }
        std::cout << ">>>>>\t" << flow << std::endl;
    }
};

template<typename cost_type>
struct OD {
    typename Graph<cost_type>::vertex_type origin;
    typename Graph<cost_type>::vertex_type destination;
    double flow;

    std::list<Path<cost_type>> paths;

    boost::numeric::ublas::compressed_matrix<double> link_flows;
    boost::numeric::ublas::compressed_matrix<double> new_link_flows;

    OD(const typename Graph<cost_type>::vertex_type& origin,
        const typename Graph<cost_type>::vertex_type& destination,
        const double& flow) : origin(origin), destination(destination), flow(flow) {}

};

template<typename cost_type>
class OD_set {
public:
    std::vector<OD<cost_type>> od_pairs;
    std::vector<std::vector<OD<cost_type>*>> ods_from_origin;
    int num_od_pairs{0};

    void add_od_pair(const typename Graph<cost_type>::vertex_type& origin,
        const typename Graph<cost_type>::vertex_type& destination,
        const double& flow, const double& budget=std::numeric_limits<double>::max()) {
        OD<cost_type> od(origin, destination, flow);
        od_pairs.push_back(od);
        num_od_pairs++;

        if (origin >= ods_from_origin.size()) {
            ods_from_origin.resize(origin + 1);
        }
        ods_from_origin[origin].push_back(&od_pairs.back());
    }
};

#endif //DEMAND_H
