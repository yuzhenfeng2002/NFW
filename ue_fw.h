//
// Created by Yuzhen Feng on 9/3/25.
//

#ifndef UE_FW_H
#define UE_FW_H

#include "graph.h"
#include "demand.h"
#include "utils.h"
#include "params.h"
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <thread>

template<typename cost_type>
class UE_FW {
public:
    Graph<cost_type>& graph;
    OD_set<cost_type>& od_set;

    UE_FW(Graph<cost_type>& graph, OD_set<cost_type>& od_set) : graph(graph), od_set(od_set) {}
    void initialization();
    void print_link_flow() { graph.print_edges(); }
    void frank_wolfe(const int& max_iter_num=100, const double& eps=1e-6);
    void print_current_iteration();
    void print_average_cost();
private:
    void update_flow_4path(Path<cost_type> path, double old_flow, double new_flow);
    double exact_line_search(const Graph<cost_type>& graph, const double& eps);
    double exact_line_search_fibonacci(const Graph<cost_type>& graph);
    double current_tstt{};
    double current_sptt{};
};

template<typename cost_type>
void FW_shortest_path(std::vector<double>& distances,
    std::vector<typename Graph<cost_type>::vertex_type>& predecessors,
    Graph<cost_type>& graph, typename Graph<cost_type>::vertex_type origin) {
    auto predecessor_map = boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index, graph.g));
    auto distance_map = boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, graph.g));
    auto weight_map = boost::get(&edge_info<cost_type>::cost, graph.g);
    boost::dijkstra_shortest_paths(graph.g, origin, boost::distance_map(distance_map).predecessor_map(predecessor_map).weight_map(weight_map));
}

template<typename cost_type>
double UE_FW<cost_type>::exact_line_search(const Graph<cost_type> &graph, const double &eps) {
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
double UE_FW<cost_type>::exact_line_search_fibonacci(const Graph<cost_type> &graph) {
    auto& n = LINE_SEARCH_F_NUM;
    auto& F = math2::F60;
    double a = 0.0, b = 1.0;
    double x1 = a + static_cast<double>(F[n - 2]) / F[n] * (b - a);
    double x2 = a + static_cast<double>(F[n - 1]) / F[n] * (b - a);

    double f1 = 0.0;
    double f2 = 0.0;

    auto edges = boost::edges(graph.g);
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
void UE_FW<cost_type>::frank_wolfe(const int& max_iter_num, const double& eps) {
    int num_iterations = 0;
    double error = std::numeric_limits<double>::max();
    std::cout << std::setw(10) << "Iteration" << std::setw(10) << "Error" << std::endl;
    int num_sources = od_set.ods_from_origin.size();
    unsigned num_threads;
    if (MULTI_THREAD) {
        num_threads = std::min(std::thread::hardware_concurrency(), static_cast<unsigned>(num_sources));
    }
    else num_threads = 1;
    size_t chunk = num_sources / num_threads;
    size_t remainder = num_sources % num_threads;
    std::vector<std::vector<Path<cost_type>>> all_origin_paths(num_sources);
    std::vector<double> all_sptt(num_sources, 0);
    double step_size = 0;
    while (num_iterations < max_iter_num && error > eps) {
        current_sptt = 0; current_tstt = 0;
        auto edges = boost::edges(graph.g);
        for (auto it = edges.first; it != edges.second; ++it) {
            auto& edge_info = graph.g[*it];
            edge_info.new_flow = 0;
            current_tstt += edge_info.flow * edge_info.cost;
        }

        std::vector<std::thread> threads;
        auto worker = [&](int start, int end) {
            for (int i = start; i < end; ++i) {
                if (od_set.ods_from_origin[i].empty()) continue;

                double sptt_of_origin = 0;
                all_origin_paths[i].clear();
                std::vector<double> distances(boost::num_vertices(graph.g));
                std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
                FW_shortest_path(distances, predecessors, graph, od_set.ods_from_origin[i][0]->origin);
                for (auto& od : od_set.ods_from_origin[i]) {
                    sptt_of_origin += od->flow * distances[od->destination];

                    Path<cost_type> path(graph);
                    path.flow = od->flow;
                    path.initialize(graph, predecessors, od->origin, od->destination);
                    path.update_cost(graph);
                    all_origin_paths[i].push_back(path);
                }
                all_sptt[i] = sptt_of_origin;
            }
        };

        size_t start_idx = 0;
        for (unsigned t = 0; t < num_threads; ++t) {
            size_t end_idx = start_idx + chunk + (t < remainder ? 1 : 0);
            threads.emplace_back(worker, start_idx, end_idx);
            start_idx = end_idx;
        }
        for (auto& t : threads) t.join();

        /* Single thread version */
        // for (int i = 0; i < num_sources; ++i) {
        //     if (od_set.ods_from_origin[i].empty()) continue;
        //
        //     double sptt_of_origin = 0;
        //     all_origin_paths[i].clear();
        //     std::vector<double> distances(boost::num_vertices(graph.g));
        //     std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
        //     FW_shortest_path(distances, predecessors, graph, od_set.ods_from_origin[i][0]->origin);
        //     for (auto& od : od_set.ods_from_origin[i]) {
        //         sptt_of_origin += od->flow * distances[od->destination];
        //
        //         Path<cost_type> path(graph);
        //         path.flow = od->flow;
        //         path.initialize(graph, predecessors, od->origin, od->destination);
        //         path.update_cost(graph);
        //         all_origin_paths[i].push_back(path);
        //     }
        //     all_sptt[i] = sptt_of_origin;
        // }

        // std::vector<double> distances(boost::num_vertices(graph.g));
        // std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
        // int last_origin_id = -1;
        //
        // for (auto& od : od_set.od_pairs) {
        //     auto origin = od.origin;
        //     auto destination = od.destination;
        //     auto flow = od.flow;
        //
        //     if (last_origin_id != origin) {
        //         FW_shortest_path(distances, predecessors, graph, origin);
        //         last_origin_id = origin;
        //     }
        //
        //     Path<cost_type> path(graph);
        //     path.flow = flow;
        //     path.initialize(graph, predecessors, origin, destination);
        //     path.update_cost(graph);
        //     current_sptt += flow * path.cost;
        //
        //     for (auto& edge : path.edge_list) {
        //         auto source = boost::source(edge, graph.g);
        //         auto target = boost::target(edge, graph.g);
        //         auto& edge_info = graph.g[edge];
        //         edge_info.new_flow += flow;
        //     }
        // }

        for (size_t i = 0; i < num_sources; ++i) {
            current_sptt += all_sptt[i];
            for (auto& path : all_origin_paths[i]) {
                for (auto& edge : path.edge_list) {
                    auto source = boost::source(edge, graph.g);
                    auto target = boost::target(edge, graph.g);
                    auto& edge_info = graph.g[edge];
                    edge_info.new_flow += path.flow;
                }
            }
        }

        error = std::abs(current_tstt - current_sptt) / current_sptt;
        num_iterations++;
        if (num_iterations % 100 == 0) {
            std::cout << std::setw(10) << num_iterations << std::setw(10) << error << std::endl;
        }

        // double step_size = 2.0 / (num_iterations + 2);
        // double step_size = exact_line_search(graph, eps * 0.001);
        step_size = exact_line_search_fibonacci(graph);

        for (auto it = edges.first; it != edges.second; ++it) {
            auto edge = *it;
            auto source = boost::source(edge, graph.g);
            auto target = boost::target(edge, graph.g);
            auto& edge_info = graph.g[edge];
            edge_info.update_flow(edge_info.flow * (1 - step_size) + edge_info.new_flow * step_size);
        }
    }
}

template<typename cost_type>
void UE_FW<cost_type>::initialization() {
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
    }
}

template<typename cost_type>
void UE_FW<cost_type>::update_flow_4path(Path<cost_type> path, double old_flow, double new_flow) {
    for (auto& edge : path.edge_list) {
        auto& edge_info = graph.g[edge];
        edge_info.update_flow(edge_info.flow - old_flow + new_flow);
    }
}

#endif //UE_FW_H
