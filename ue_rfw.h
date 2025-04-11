//
// Created by Yuzhen Feng on 3/4/25.
//

#ifndef UE_RFW_H
#define UE_RFW_H

#include "graph.h"
#include "demand.h"
#include "mqcf.h"
#include "cost.h"
#include "utils.h"
#include "params.h"
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <thread>

template<typename cost_type>
void RFW_shortest_path(std::vector<double>& distances,
    std::vector<typename Graph<cost_type>::vertex_type>& predecessors,
    Graph<cost_type>& graph, typename Graph<cost_type>::vertex_type origin) {
    auto predecessor_map = boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index, graph.g));
    auto distance_map = boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, graph.g));
    auto weight_map = boost::get(&edge_info<cost_type>::cost, graph.g);
    boost::dijkstra_shortest_paths(graph.g, origin, boost::distance_map(distance_map).predecessor_map(predecessor_map).weight_map(weight_map));
}

template<typename cost_type>
class UE_RFW {
public:
    Graph<cost_type>& graph;
    OD_set<cost_type>& od_set;
    std::pair<typename boost::graph_traits<typename Graph<cost_type>::graph_type>::edge_iterator,
                typename boost::graph_traits<typename Graph<cost_type>::graph_type>::edge_iterator> edges;

    UE_RFW(Graph<cost_type>& graph, OD_set<cost_type>& od_set) : graph(graph), od_set(od_set) {edges = boost::edges(graph.g);}
    void initialization();
    void print_link_flow() { graph.print_edges(); }
    void revised_frank_wolfe(const int& max_iter_num=100, const double& eps=1e-6);
private:
    void update_flow_4path(Path<cost_type> path, double old_flow, double new_flow);
    double exact_line_search_fibonacci(std::vector<double> alphas);
    double current_tstt{};
    double current_sptt{};
};

template<typename cost_type>
double UE_RFW<cost_type>::exact_line_search_fibonacci(std::vector<double> alphas) {
    auto& n = LINE_SEARCH_F_NUM;
    auto& F = math2::F60;
    double a = 0.0, b = 1.0;
    double x1 = a + static_cast<double>(F[n - 2]) / F[n] * (b - a);
    double x2 = a + static_cast<double>(F[n - 1]) / F[n] * (b - a);

    double f1 = 0.0;
    double f2 = 0.0;

    for (auto it = edges.first; it != edges.second; ++it) {
        auto& edge_info = graph.g[*it];
        double new_flow_mid1 = 0;
        double new_flow_mid2 = 0;
        for (int i = 0; i < alphas.size(); ++i) {
            new_flow_mid1 += edge_info.orgin_flows[i] * (1 - x1 * alphas[i]) + edge_info.origin_new_flows[i] * x1 * alphas[i];
        }
        for (int i = 0; i < alphas.size(); ++i) {
            new_flow_mid2 += edge_info.orgin_flows[i] * (1 - x2 * alphas[i]) + edge_info.origin_new_flows[i] * x2 * alphas[i];
        }
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
                double new_flow_mid2 = 0;
                for (int i = 0; i < alphas.size(); ++i) {
                    new_flow_mid2 += edge_info.orgin_flows[i] * (1 - x2 * alphas[i]) + edge_info.origin_new_flows[i] * x2 * alphas[i];
                }
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
                double new_flow_mid1 = 0;
                for (int i = 0; i < alphas.size(); ++i) {
                    new_flow_mid1 += edge_info.orgin_flows[i] * (1 - x1 * alphas[i]) + edge_info.origin_new_flows[i] * x1 * alphas[i];
                }
                f1 += edge_info.cost_fun.integral(new_flow_mid1);
            }
        }
    }

    return (a + b) / 2.0;
}

template<typename cost_type>
void UE_RFW<cost_type>::revised_frank_wolfe(const int& max_iter_num, const double& eps) {
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
    std::vector<double> all_sptt(num_sources, 0);
    std::vector<double> all_tstt(num_sources, 0);
    double step_size = 0;
    std::vector<double> alphas(num_sources, 0);
    while (num_iterations < max_iter_num && error > eps) {
        current_sptt = 0; current_tstt = 0;

        for (auto it = edges.first; it != edges.second; ++it) {
            auto edge = *it;
            auto source = boost::source(edge, graph.g);
            auto target = boost::target(edge, graph.g);
            auto& edge_info = graph.g[edge];
            double new_flow = 0;
            for (int i = 0; i < od_set.ods_from_origin.size(); ++i) {
                new_flow += edge_info.orgin_flows[i] * (1 - step_size * alphas[i]) + edge_info.origin_new_flows[i] * step_size * alphas[i];
                edge_info.orgin_flows[i] = edge_info.orgin_flows[i] * (1 - step_size * alphas[i]) + edge_info.origin_new_flows[i] * step_size * alphas[i];
                edge_info.origin_new_flows[i] = 0;
            }
            edge_info.update_flow(new_flow);
            edge_info.new_flow = 0;
            current_tstt += edge_info.flow * edge_info.cost;
        }

        std::vector<std::thread> threads;
        auto worker = [&](int start, int end) {
            for (int i = start; i < end; ++i) {
                if (od_set.ods_from_origin[i].empty()) continue;

                double sptt_of_origin = 0;
                double tstt_of_origin = 0;
                double denominator = 0;
                std::vector<double> distances(boost::num_vertices(graph.g));
                std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
                RFW_shortest_path(distances, predecessors, graph, od_set.ods_from_origin[i][0]->origin);
                for (auto& od : od_set.ods_from_origin[i]) {
                    Path<cost_type> path(graph);
                    path.initialize(graph, predecessors, od->origin, od->destination);
                    for (auto& edge : path.edge_list) {
                        auto source = boost::source(edge, graph.g);
                        auto target = boost::target(edge, graph.g);
                        auto& edge_info = graph.g[edge];
                        edge_info.origin_new_flows[i] += od->flow;
                    }
                }
                for (auto it = edges.first; it != edges.second; ++it) {
                    auto edge = *it;
                    auto source = boost::source(edge, graph.g);
                    auto target = boost::target(edge, graph.g);
                    auto& edge_info = graph.g[edge];
                    tstt_of_origin += edge_info.orgin_flows[i] * edge_info.cost;
                    sptt_of_origin += edge_info.origin_new_flows[i] * edge_info.cost;
                    denominator += edge_info.cost_derivative * ((edge_info.orgin_flows[i] - edge_info.origin_new_flows[i])
                    * (edge_info.orgin_flows[i] - edge_info.origin_new_flows[i]));
                }
                all_sptt[i] = sptt_of_origin;
                all_tstt[i] = tstt_of_origin;
                if (denominator < 1e-20) {
                    alphas[i] = 1;
                }
                else {
                    alphas[i] = (tstt_of_origin - sptt_of_origin) / denominator;
                    if (alphas[i] > 1) {
                        alphas[i] = 1;
                    }
                }
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
        //     double tstt_of_origin = 0;
        //     double denominator = 0;
        //     std::vector<double> distances(boost::num_vertices(graph.g));
        //     std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));
        //     RFW_shortest_path(distances, predecessors, graph, od_set.ods_from_origin[i][0]->origin);
        //     for (auto& od : od_set.ods_from_origin[i]) {
        //         // sptt_of_origin += od->flow * distances[od->destination];
        //         Path<cost_type> path(graph);
        //         path.initialize(graph, predecessors, od->origin, od->destination);
        //         for (auto& edge : path.edge_list) {
        //             auto source = boost::source(edge, graph.g);
        //             auto target = boost::target(edge, graph.g);
        //             auto& edge_info = graph.g[edge];
        //             edge_info.origin_new_flows[i] += od->flow;
        //         }
        //     }
        //     for (auto it = edges.first; it != edges.second; ++it) {
        //         auto edge = *it;
        //         auto source = boost::source(edge, graph.g);
        //         auto target = boost::target(edge, graph.g);
        //         auto& edge_info = graph.g[edge];
        //         tstt_of_origin += edge_info.orgin_flows[i] * edge_info.cost;
        //         sptt_of_origin += edge_info.origin_new_flows[i] * edge_info.cost;
        //         denominator += edge_info.cost_derivative * ((edge_info.orgin_flows[i] - edge_info.origin_new_flows[i])
        //         * (edge_info.orgin_flows[i] - edge_info.origin_new_flows[i]));
        //     }
        //     all_sptt[i] = sptt_of_origin;
        //     all_tstt[i] = tstt_of_origin;
        //     if (denominator < 1e-20) {
        //         alphas[i] = 1;
        //     }
        //     else {
        //         alphas[i] = (tstt_of_origin - sptt_of_origin) / denominator;
        //         if (alphas[i] > 1) {
        //             alphas[i] = 1;
        //         }
        //     }
        // }

        for (size_t i = 0; i < num_sources; ++i) {
            current_sptt += all_sptt[i];
        }

        error = std::abs(current_tstt - current_sptt) / current_sptt;
        num_iterations++;
        if (num_iterations % 100 == 0) {
            std::cout << std::setw(10) << num_iterations << std::setw(10) << error << std::endl;
        }

        // step_size = 2.0 / (num_iterations + 2);
        step_size = exact_line_search_fibonacci(alphas);
    }
}

template<typename cost_type>
void UE_RFW<cost_type>::initialization() {
    std::vector<double> distances(boost::num_vertices(graph.g));
    std::vector<typename Graph<cost_type>::vertex_type> predecessors(boost::num_vertices(graph.g));

    int num_of_sources = od_set.ods_from_origin.size();
    for (auto it = edges.first; it != edges.second; ++it) {
        auto edge = *it;
        auto& edge_info = graph.g[edge];
        edge_info.origin_new_flows.resize(num_of_sources);
        edge_info.orgin_flows.resize(num_of_sources);
    }

    for (size_t i = 0; i < od_set.ods_from_origin.size(); ++i) {
        if (od_set.ods_from_origin[i].empty()) continue;
        RFW_shortest_path(distances, predecessors, graph, od_set.ods_from_origin[i][0]->origin);
        for (auto& od : od_set.ods_from_origin[i]) {
            Path<cost_type> path(graph);
            path.flow = od->flow;
            path.initialize(graph, predecessors, od->origin, od->destination);
            update_flow_4path(path, 0, od->flow);
            for (auto& edge : path.edge_list) {
                auto& edge_info = graph.g[edge];
                edge_info.orgin_flows[i] += od->flow;
            }
        }
    }
}

template<typename cost_type>
void UE_RFW<cost_type>::update_flow_4path(Path<cost_type> path, double old_flow, double new_flow) {
    for (auto& edge : path.edge_list) {
        auto& edge_info = graph.g[edge];
        edge_info.update_flow(edge_info.flow - old_flow + new_flow);
    }
}

#endif //UE_RFW_H
