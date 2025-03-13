//
// Created by Yuzh on 9/3/25.
//

#ifndef UE_NFW_H
#define UE_NFW_H

#include "graph.h"
#include "demand.h"
#include "mqcf.h"
#include "cost.h"
#include "utils.h"
#include "params.h"
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
    std::pair<typename boost::graph_traits<typename Graph<cost_type>::graph_type>::edge_iterator,
                typename boost::graph_traits<typename Graph<cost_type>::graph_type>::edge_iterator> edges;

    UE_NFW(Graph<cost_type>& graph, OD_set<cost_type>& od_set) : graph(graph), od_set(od_set) {edges = boost::edges(graph.g);}
    void initialization();
    void print_link_flow() { graph.print_edges(); }
    void newton_frank_wolfe(const int& max_iter_num=100, const double& eps=1e-6);
private:
    void update_flow_4path(Path<cost_type> path, double old_flow, double new_flow);
    double exact_line_search(const Graph<cost_type>& graph, const double& eps);
    double exact_line_search_fibonacci(const Graph<cost_type>& graph, const double& eps);
    double current_tstt{};
    double current_sptt{};
};

template<typename cost_type>
double UE_NFW<cost_type>::exact_line_search(const Graph<cost_type> &graph, const double &eps) {
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
double UE_NFW<cost_type>::exact_line_search_fibonacci(const Graph<cost_type> &graph, const double &eps) {
    auto& n = LINE_SEARCH_F_NUM;
    auto& F = math2::F60;
    double a = 0.0, b = 1.0;
    double x1 = a + static_cast<double>(F[n - 2]) / F[n] * (b - a);
    double x2 = a + static_cast<double>(F[n - 1]) / F[n] * (b - a);

    double f1 = 0.0;
    double f2 = 0.0;

    for (auto it = edges.first; it != edges.second; ++it) {
        auto& edge_info = graph.g[*it];
        double new_flow_mid1 = edge_info.flow * (1 - x1) + edge_info.new_flow * x1;
        double new_flow_mid2 = edge_info.flow * (1 - x2) + edge_info.new_flow * x2;
        f1 += edge_info.cost_fun.integral(new_flow_mid1);
        f2 += edge_info.cost_fun.integral(new_flow_mid2);
    }

    for (int k = n; k > 1; --k) {
        if (f1 > f2) { // Minimum is in [x1, b]
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + static_cast<double>(F[k - 1]) / F[k] * (b - a);
            f2 = 0;
            for (auto it = edges.first; it != edges.second; ++it) {
                auto& edge_info = graph.g[*it];
                double new_flow_mid2 = edge_info.flow * (1 - x2) + edge_info.new_flow * x2;
                f2 += edge_info.cost_fun.integral(new_flow_mid2);
            }
        } else { // Minimum is in [a, x2]
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + static_cast<double>(F[k - 2]) / F[k] * (b - a);
            f1 = 0;
            for (auto it = edges.first; it != edges.second; ++it) {
                auto& edge_info = graph.g[*it];
                double new_flow_mid1 = edge_info.flow * (1 - x1) + edge_info.new_flow * x1;
                f1 += edge_info.cost_fun.integral(new_flow_mid1);
            }
        }
    }

    return (a + b) / 2.0;
}

template<typename cost_type>
void UE_NFW<cost_type>::newton_frank_wolfe(const int& max_iter_num, const double& eps) {
    int num_iterations = 0;
    double error = std::numeric_limits<double>::max();
    std::vector<double> distances(boost::num_vertices(graph.g));
    std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
    std::cout << std::setw(10) << "Iteration" << std::setw(10) << "Error" << std::endl;
    while (num_iterations < max_iter_num && error > eps) {
        current_sptt = 0; current_tstt = 0;

        for (int i = 0; i < od_set.ods_from_origin.size(); ++i) {
            if (od_set.ods_from_origin[i].empty()) continue;
            MQCF<cost_type> mqcf(graph, od_set, i);

            /* For debugging start*/
            // Graph<quad> subgraph;
            // subgraph.num_edges = boost::num_vertices(graph.g);
            // subgraph.num_vertices = boost::num_vertices(graph.g);
            // auto vertices = boost::vertices(graph.g);
            // for (auto it = vertices.first; it != vertices.second; ++it) {
            //     auto vertex = *it;
            //     boost::add_vertex(subgraph.g);
            // }
            // auto edges = boost::edges(graph.g);
            // for (auto it = edges.first; it != edges.second; ++it) {
            //     auto edge = *it;
            //     auto source = boost::source(edge, graph.g);
            //     auto target = boost::target(edge, graph.g);
            //     boost::add_edge(source, target, subgraph.g);
            // }
            // auto sub_edges = boost::edges(subgraph.g);
            // for (auto it = sub_edges.first; it != sub_edges.second; ++it) {
            //     auto edge = *it;
            //     auto source = boost::source(edge, graph.g);
            //     auto target = boost::target(edge, graph.g);
            //     auto& edge_info = graph.g[edge];
            //     auto sub_edge = boost::edge(source, target, subgraph.g).first;
            //     auto& sub_edge_info = subgraph.g[sub_edge];
            //     auto sub_mqcf_edge = boost::edge(source, target, mqcf.QCgraph).first;
            //     auto& sub_mqcf_edge_info = mqcf.QCgraph[sub_mqcf_edge];
            //     sub_edge_info.cost_fun.initialize(sub_mqcf_edge_info.c, sub_mqcf_edge_info.d);
            // }
            // OD_set<quad> sub_od_set;
            // for (auto& od : od_set.ods_from_origin[i]) {
            //     sub_od_set.add_od_pair(od->origin, od->destination, od->flow);
            // }
            // UE_GP<quad> sub_ue_fw(subgraph, sub_od_set);
            // sub_ue_fw.initialization();
            // sub_ue_fw.gradient_projection(10000, eps * 0.1);
            // sub_ue_fw.print_link_flow();
            /* For debugging end*/

            mqcf.basic_algorithm(1000000, eps * 0.0001);
            auto qc_edges = boost::edges(mqcf.QCgraph);

            for (auto it = qc_edges.first; it != qc_edges.second; ++it) {
                auto& qc_edge_info = mqcf.QCgraph[*it];
                auto source = boost::source(*it, mqcf.QCgraph);
                auto target = boost::target(*it, mqcf.QCgraph);

                auto edge = boost::edge(source, target, graph.g).first;
                auto& edge_info = graph.g[edge];
                edge_info.new_flow += qc_edge_info.flow;
                edge_info.origin_new_flows[i] = qc_edge_info.flow;
            }
        }

        // double step_size = 2.0 / (num_iterations + 2);
        double step_size = exact_line_search_fibonacci(graph, eps * 0.0001);

        for (auto it = edges.first; it != edges.second; ++it) {
            auto edge = *it;
            auto source = boost::source(edge, graph.g);
            auto target = boost::target(edge, graph.g);
            auto& edge_info = graph.g[edge];
            edge_info.update_flow(edge_info.flow * (1 - step_size) + edge_info.new_flow * step_size);
            for (int i = 0; i < od_set.ods_from_origin.size(); ++i) {
                edge_info.orgin_flows[i] = edge_info.orgin_flows[i] * (1 - step_size) + edge_info.origin_new_flows[i] * step_size;
                edge_info.origin_new_flows[i] = 0;
            }
            current_tstt += edge_info.flow * edge_info.cost;
            edge_info.new_flow = 0;
        }

        for (size_t i = 0; i < od_set.ods_from_origin.size(); ++i) {
            if (od_set.ods_from_origin[i].empty()) continue;
            NFW_shortest_path(distances, predecessors, graph, i);
            for (auto& od : od_set.ods_from_origin[i]) {
                current_sptt += od->flow * distances[od->destination];
            }
        }

        error = std::abs(current_tstt - current_sptt) / current_sptt;
        num_iterations++;
        if (num_iterations % 1 == 0) {
            std::cout << std::setw(10) << num_iterations << std::setw(10) << error << std::endl;
        }
    }
}

template<typename cost_type>
void UE_NFW<cost_type>::initialization() {
    std::vector<double> distances(boost::num_vertices(graph.g));
    std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));

    for (auto it = edges.first; it != edges.second; ++it) {
        auto edge = *it;
        auto& edge_info = graph.g[edge];
        edge_info.origin_new_flows.resize(graph.num_vertices);
        edge_info.orgin_flows.resize(graph.num_vertices);
    }

    for (size_t i = 0; i < od_set.ods_from_origin.size(); ++i) {
        if (od_set.ods_from_origin[i].empty()) continue;
        NFW_shortest_path(distances, predecessors, graph, i);
        for (auto& od : od_set.ods_from_origin[i]) {
            Path<cost_type> path(graph);
            path.flow = od->flow;
            path.initialize(graph, predecessors, i, od->destination);
            update_flow_4path(path, 0, od->flow);
            for (auto& edge : path.edge_list) {
                auto& edge_info = graph.g[edge];
                edge_info.orgin_flows[i] += od->flow;
            }
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
