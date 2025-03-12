//
// Created by Yuzhen Feng on 19/2/2025.
//

#ifndef UE_GP_H
#define UE_GP_H

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <iostream>
#include "graph.h"
#include "demand.h"

template<typename cost_type>
class UE_GP {
public:
    Graph<cost_type>& graph;
    OD_set<cost_type>& od_set;

    UE_GP(Graph<cost_type>& graph, OD_set<cost_type>& od_set) : graph(graph), od_set(od_set) {}
    void initialization();
    void print_link_flow() { graph.print_edges(); }
    void gradient_projection(const int& max_iter_num=100, const double& eps=1e-6, const double& step_size=0.25);
    void print_current_iteration();
    void print_average_cost();
private:
    void update_flow_4path(Path<cost_type> path, double old_flow, double new_flow);
    double current_tstt{};
    double current_sptt{};
};

template<typename cost_type>
void UE_GP<cost_type>::update_flow_4path(Path<cost_type> path, double old_flow, double new_flow) {
    for (auto& edge : path.edge_list) {
        auto edge_info = &graph.g[edge];
        edge_info->update_flow(edge_info->flow - old_flow + new_flow);
    }
}

template<typename cost_type>
void UE_GP<cost_type>::gradient_projection(const int& max_iter_num, const double& eps, const double& step_size) {
    int num_iterations = 0;
    double error = std::numeric_limits<double>::max();
    std::cout << std::setw(10) << "Iteration" << std::setw(10) << "Error" << std::endl;
    while (num_iterations < max_iter_num && error > eps) {
        std::vector<double> distances(graph.num_vertices);
        std::vector<typename Graph<cost_type>::vertex_type> predecessors(graph.num_vertices);
        current_sptt = 0; current_tstt = 0;
        for (auto& od : od_set.od_pairs) {
            shortest_path(distances, predecessors, graph, od.origin);

            Path<cost_type> new_path(graph);
            new_path.initialize(graph, predecessors, od.origin, od.destination);
            new_path.cost = distances[od.destination];

            bool is_new_path_generated = true;
            auto old_same_path_iter = od.paths.end();
            double current_total_flow = 0;

            current_sptt += od.flow * new_path.cost;

            for (auto path_it = od.paths.begin(); path_it != od.paths.end(); ) {
                auto& path = *path_it;
                path.update_cost(graph);
                current_tstt += path.flow * path.cost;
                ++path_it;
            }

            for (auto path_it = od.paths.begin(); path_it != od.paths.end(); ) {
                auto& path = *path_it;
                double first_order_derivative = path.cost - new_path.cost;
                double second_order_derivative = 0;

                bool is_same_path = true;
                auto edges = boost::edges(graph.g);
                for (auto it = edges.first; it != edges.second; ++it) {
                    auto edge = *it;
                    auto head = boost::source(edge, graph.g);
                    auto tail = boost::target(edge, graph.g);

                    bool is_same_link = path.edge_matrix(head, tail) == new_path.edge_matrix(head, tail);
                    is_same_path = is_same_path && is_same_link;
                    second_order_derivative += graph.g[edge].cost_derivative * static_cast<int>(!is_same_link);
                }

                is_new_path_generated = is_new_path_generated && !is_same_path;
                if (!is_same_path) {
                    double path_old_flow = path.flow;
                    path.flow -= step_size * first_order_derivative / second_order_derivative;
                    if (path.flow < 0) {
                        path.flow = 0;
                        update_flow_4path(path, path_old_flow, 0);
                        auto next_path_it = std::next(path_it);
                        od.paths.erase(path_it);
                        path_it = next_path_it;
                        continue;
                    }
                    update_flow_4path(path, path_old_flow, path.flow);
                    current_total_flow += path.flow;
                }
                else {
                    old_same_path_iter = path_it;
                }
                ++path_it;
            }

            if (is_new_path_generated) {
                new_path.flow = od.flow - current_total_flow;
                od.paths.push_back(new_path);
                update_flow_4path(new_path, 0, new_path.flow);
            }
            else {
                double new_path_old_flow = old_same_path_iter->flow;
                old_same_path_iter->flow = od.flow - current_total_flow;
                update_flow_4path(*old_same_path_iter, new_path_old_flow, old_same_path_iter->flow);
            }
        }
        num_iterations++;
        error = std::abs(current_tstt - current_sptt) / current_sptt;
        if (num_iterations % 1 == 0) {
            std::cout << std::setw(10) << num_iterations << std::setw(10) << error << std::endl;
        }
    }
}

template<typename cost_type>
void shortest_path(std::vector<double>& distances,
    std::vector<typename Graph<cost_type>::vertex_type>& predecessors,
    Graph<cost_type>& graph, typename Graph<cost_type>::vertex_type origin) {
    auto predecessor_map = boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index, graph.g));
    auto distance_map = boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, graph.g));
    auto weight_map = boost::get(&edge_info<cost_type>::cost, graph.g);
    boost::dijkstra_shortest_paths(graph.g, origin, boost::distance_map(distance_map).predecessor_map(predecessor_map).weight_map(weight_map));

    /* Bellman-Ford algorithm for debugging start */
    // for (size_t i = 0; i < distances.size(); ++i) {
    //     distances[i] = std::numeric_limits<double>::max();
    // }
    // distances[origin] = 0;
    // auto vertices = boost::vertices(graph.g);
    // auto edges = boost::edges(graph.g);
    // int n = boost::num_vertices(graph.g);
    // for (size_t i = 0; i < n; ++i) {
    //     for (auto eit = edges.first; eit != edges.second; ++eit) {
    //         auto e = *eit;
    //         auto& edge_info = graph.g[e];
    //         auto u = source(e, graph.g);
    //         auto v = target(e, graph.g);
    //         double cost = edge_info.cost;
    //         if (distances[u] + cost < distances[v]) {
    //             distances[v] = distances[u] + cost;
    //             predecessors[v] = u;
    //         }
    //     }
    // }
    /* Bellman-Ford algorithm for debugging end */
}

template<typename cost_type>
void UE_GP<cost_type>::initialization() {
    std::vector<double> distances(boost::num_vertices(graph.g));
    std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
    int last_origin_id = -1;
    for (auto& od : od_set.od_pairs) {
        auto origin = od.origin;
        auto destination = od.destination;
        auto flow = od.flow;
        if (last_origin_id != origin) {
            shortest_path(distances, predecessors, graph, origin);
            last_origin_id = origin;
        }
        Path<cost_type> path(graph);
        path.flow = flow;
        path.initialize(graph, predecessors, origin, destination);
        od.paths.push_back(path);
        update_flow_4path(path, 0, flow);
    }
}

template<typename cost_type>
void UE_GP<cost_type>::print_current_iteration() {
    for (auto& od : od_set.od_pairs) {
        for (auto& path : od.paths) {
            path.print_path();
        }
    }
    graph.print_edges();
}

template<typename cost_type>
void UE_GP<cost_type>::print_average_cost() {
    std::vector<std::vector<double>> avg_cost(graph.num_vertices, std::vector<double>(graph.num_vertices, 0.));
    for (auto& od : od_set.od_pairs) {
        for (auto& path : od.paths) {
            path.update_cost(graph);
            avg_cost[od.origin][od.destination] += path.cost * path.flow;
        }
        if (od.paths.empty()) {
            continue;
        }
        avg_cost[od.origin][od.destination] /= od.flow;
    }

    for (size_t i = 0; i < avg_cost.size(); ++i) {
        std::cout << "Origin  " << i + 1 << std::endl;
        int count = 0;
        for (size_t j = 0; j < avg_cost[i].size(); ++j) {
            std::cout << std::setw(5) << std::right << j + 1 << " : "
            << std::fixed << std::setprecision(1) << std::setw(8) << std::right << avg_cost[i][j] << "; ";
            count++;
            if (count % 5 == 0) {
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;
    }
}

#endif //UE_GP_H
