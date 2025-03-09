//
// Created by Yuzh on 9/3/25.
//

#ifndef UE_NFW_H
#define UE_NFW_H

#include "graph.h"
#include "demand.h"
#include "mqcf.h"
#include <boost/graph/dijkstra_shortest_paths.hpp>

template<typename cost_type>
void NFW_shortest_path(std::vector<double>& distances,
    std::vector<typename Graph<cost_type>::vertex_type>& predecessors,
    Graph<cost_type>& graph, typename Graph<cost_type>::vertex_type origin) {
    auto predecessor_map = boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index, graph.g));
    auto distance_map = boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, graph.g));
    auto weight_map = boost::get(&edge_info<cost_type>::cost, graph.g);
    boost::dijkstra_shortest_paths(graph.g, origin, boost::distance_map(distance_map).predecessor_map(predecessor_map).weight_map(weight_map));
}

template<typename cost_type>
class UE_NFW {
public:
    Graph<cost_type>& graph;
    OD_set<cost_type>& od_set;

    UE_NFW(Graph<cost_type>& graph, OD_set<cost_type>& od_set) : graph(graph), od_set(od_set) {}
    void initialization();
    void print_link_flow() { graph.print_edges(); }
    void newton_frank_wolfe(const int& max_iter_num=100, const double& eps=1e-6);
private:
    void update_flow_4path(Path<cost_type> path, double old_flow, double new_flow);
    double exact_line_search(const Graph<cost_type>& graph, const double& eps);
    double current_tstt{};
    double current_sptt{};
};



template<typename cost_type>
double UE_NFW<cost_type>::exact_line_search(const Graph<cost_type> &graph, const double &eps) {
    auto edges = boost::edges(graph.g);
    double low = 0.0;
    double high = 1.0;
    double mid = 0.0;

    while ((high - low) > eps) {
        mid = (low + high) / 2.0;
        double total_cost_integral_low = 0.0;
        double total_cost_integral_high = 0.0;

        for (auto it = edges.first; it != edges.second; ++it) {
            auto edge = *it;
            auto& edge_info = graph.g[edge];

            double new_flow_low = edge_info.flow * (1 - low) + edge_info.new_flow * low;
            double new_flow_high = edge_info.flow * (1 - high) + edge_info.new_flow * high;

            total_cost_integral_low += edge_info.cost_fun.integral(new_flow_low);
            total_cost_integral_high += edge_info.cost_fun.integral(new_flow_high);
        }

        if (total_cost_integral_low < total_cost_integral_high) {
            high = mid;
        } else {
            low = mid;
        }
    }

    return mid;
}


template<typename cost_type>
void UE_NFW<cost_type>::newton_frank_wolfe(const int& max_iter_num, const double& eps) {
    int num_iterations = 0;
    double error = std::numeric_limits<double>::max();
    std::cout << std::setw(10) << "Iteration" << std::setw(10) << "Error" << std::endl;
    while (num_iterations < max_iter_num && error > eps) {
        auto edges = boost::edges(graph.g);
        for (auto it = edges.first; it != edges.second; ++it) {
            auto& edge_info = graph.g[*it];
            edge_info.new_flow = 0;
        }

        std::vector<double> distances(boost::num_vertices(graph.g));
        std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
        int last_origin_id = -1;

        current_sptt = 0; current_tstt = 0;
        for (auto& od : od_set.od_pairs) {
            auto origin = od.origin;
            auto destination = od.destination;
            auto flow = od.flow;

            if (last_origin_id != origin) {
                NFW_shortest_path(distances, predecessors, graph, origin);
                last_origin_id = origin;
            }

            Path<cost_type> path(graph);
            path.flow = flow;
            path.initialize(graph, predecessors, origin, destination);
            path.update_cost(graph);
            current_sptt += flow * path.cost;
        }

        print_link_flow();
        for (int i = 0; i < od_set.ods_from_origin.size(); ++i) {
            MQCF<cost_type> mqcf(graph, od_set, i);
            mqcf.basic_algorithm(100, 1e-6);
            auto qc_edges = boost::edges(mqcf.QCgraph);
            for (auto it = qc_edges.first; it != qc_edges.second; ++it) {
                auto& qc_edge_info = mqcf.QCgraph[*it];
                auto source = boost::source(*it, mqcf.QCgraph);
                auto target = boost::target(*it, mqcf.QCgraph);
                if (!qc_edge_info.is_reverse) {
                    auto edge = boost::edge(source, target, graph.g).first;
                    auto& edge_info = graph.g[edge];
                    edge_info.new_flow = qc_edge_info.flow;
                }
            }
        }

        // double step_size = 2.0 / (num_iterations + 2);
        double step_size = exact_line_search(graph, eps);

        for (auto it = edges.first; it != edges.second; ++it) {
            auto edge = *it;
            auto source = boost::source(edge, graph.g);
            auto target = boost::target(edge, graph.g);
            auto& edge_info = graph.g[edge];
            edge_info.update_flow(edge_info.flow * (1 - step_size) + edge_info.new_flow * step_size);
            current_tstt += edge_info.flow * edge_info.cost;
        }

        for (auto& od : od_set.od_pairs) {
            od.link_flows = od.link_flows * (1 - step_size) + od.new_link_flows * step_size;
        }

        error = std::abs(current_tstt - current_sptt) / current_sptt;
        num_iterations++;
        if (num_iterations % 100 == 0) {
            std::cout << std::setw(10) << num_iterations << std::setw(10) << error << std::endl;
        }
    }
}

template<typename cost_type>
void UE_NFW<cost_type>::initialization() {
    std::vector<double> distances(boost::num_vertices(graph.g));
    std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
    int last_origin_id = -1;

    for (auto& od : od_set.od_pairs) {
        auto origin = od.origin;
        auto destination = od.destination;
        auto flow = od.flow;
        if (last_origin_id != origin) {
            FW_shortest_path(distances, predecessors, graph, origin);
            last_origin_id = origin;
        }
        Path<cost_type> path(graph);
        path.flow = flow;
        path.initialize(graph, predecessors, origin, destination);
        update_flow_4path(path, 0, flow);

        od.link_flows.resize(graph.num_vertices, graph.num_vertices, false);
        od.new_link_flows.resize(graph.num_vertices, graph.num_vertices, false);

        for (auto& edge : path.edge_list) {
            auto source = boost::source(edge, graph.g);
            auto target = boost::target(edge, graph.g);
            od.link_flows(source, target) += flow;
        }
    }
}

template<typename cost_type>
void UE_NFW<cost_type>::update_flow_4path(Path<cost_type> path, double old_flow, double new_flow) {
    for (auto& edge : path.edge_list) {
        auto& edge_info = graph.g[edge];
        edge_info.update_flow(edge_info.flow - old_flow + new_flow);
    }
}

#endif //UE_NFW_H
