//
// Created by Yuzhen Feng on 19/2/2025.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <iomanip>

// Define the vertex information structure
struct vertex_info {
    bool centroid = true;
};

// Define the edge information structure template
template<typename cost_type> struct edge_info {
    double flow;
    double cost;
    double cost_derivative;
    cost_type cost_fun;

    double new_flow{0.};

    std::vector<double> orgin_flows;
    std::vector<double> origin_new_flows;

    edge_info() : flow(0.), cost(0.), cost_derivative(0.) {}

    void update_only_flow(const double& new_flow) {
        flow = new_flow;
    }

    void update_only_fun() {
        this->cost_fun.update(flow, this->cost, this->cost_derivative);
    }

    void update_flow(const double& new_flow) {
        update_only_flow(new_flow);
        if (flow < 0) flow = 0;
        update_only_fun();
    }
};

// Define the Graph class template
template<typename cost_type>
class Graph {
public:
    typedef boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS, vertex_info,
    boost::property<boost::edge_index_t, std::size_t, edge_info<cost_type>>> graph_type;
    graph_type g;
    typedef typename graph_type::vertex_descriptor vertex_type;
    typedef typename graph_type::edge_descriptor edge_type;
    int num_vertices{};
    int num_edges{};
    void print_vertices();
    void print_edges();
};

template<typename cost_type>
void Graph<cost_type>::print_vertices() {
    std::cout << std::left << std::setw(10) << "Vertices: " << num_vertices << std::endl;
    std::cout << std::left << std::setw(10) << "ID" << std::left << std::setw(10) << "Centroid" << std::endl;
    auto vertices = boost::vertices(g);
    for (auto it = vertices.first; it != vertices.second; ++it) {
        std::cout << std::left << std::setw(10) << *it << std::left << std::setw(10) << g[*it].centroid << std::endl;
    }
}

template<typename cost_type>
void Graph<cost_type>::print_edges() {
    std::cout << "Edges: " << num_edges << std::endl;
    std::cout << std::left << std::setw(10) << "Source"
              << std::left << std::setw(10) << "Target"
              << std::left << std::setw(10) << "Flow"
              << std::left << std::setw(10) << "Cost"
              << std::left << std::setw(20) << "Cost Derivative" << std::endl;
    auto edges = boost::edges(g);
    for (auto it = edges.first; it != edges.second; ++it) {
        auto source = boost::source(*it, g);
        auto target = boost::target(*it, g);
        const auto& edge_info = g[*it];
        std::cout << std::left << std::setw(10) << source
                  << std::left << std::setw(10) << target
                  << std::left << std::setw(10) << edge_info.flow
                  << std::left << std::setw(10) << edge_info.cost
                  << std::left << std::setw(20) << edge_info.cost_derivative << std::endl;
    }
}

#endif //GRAPH_H
