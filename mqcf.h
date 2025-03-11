//
// Created by 封钰震 on 9/3/25.
//

#ifndef MQCF_H
#define MQCF_H

#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <cmath>
#include "graph.h"

struct NodeProps {
    double demand{0};      // Node demand (b_i)
    double excess{0};      // Current excess (ρ_f(i) - b_i)
    double pi{0};          // Dual variable (π_i)
};

struct EdgeProps {
    double c{0};           // Quadratic coefficient C_ij(f) = c*f² + d*f
    double d{0};           // Linear coefficient
    double flow{0};        // Current flow f_ij
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, NodeProps, EdgeProps> QCGraph;
typedef boost::graph_traits<QCGraph>::vertex_descriptor QCVertex;
typedef boost::graph_traits<QCGraph>::edge_descriptor QCEdge;

template<typename cost_type>
class MQCF {
public:
    QCGraph QCgraph;
    unsigned long origin;
    std::vector<std::pair<int, double>> destinations;

    MQCF(const Graph<cost_type>& graph, const OD_set<cost_type>& od_set, int& origin);

    void basic_algorithm(int max_iter=100, double epsilon=1e-6);
private:
    double compute_initial_delta();
    void adjust_flow(const double& delta);
    void augment_flow(const std::vector<std::pair<QCEdge, bool>>& pred, double delta, QCVertex& s, QCVertex& t);
    void update_potentials(const std::vector<double>& dist);
    void print_edges();
    void print_vertices();
};

template<typename cost_type>
MQCF<cost_type>::MQCF(const Graph<cost_type>& graph, const OD_set<cost_type>& od_set, int& origin) {
    this->origin = origin;

    for (auto vp = boost::vertices(graph.g); vp.first != vp.second; ++vp.first) {
        auto v = *vp.first;
        boost::add_vertex(QCgraph);
    }

    for (auto od : od_set.ods_from_origin[origin]) {
        auto destination = od->destination;
        auto flow = od->flow;
        destinations.push_back(std::make_pair(destination, flow));
        QCgraph[destination].demand += flow;
        QCgraph[destination].excess -= flow;
        QCgraph[origin].demand -= flow;
        QCgraph[origin].excess += flow;
    }

    for (auto ep = boost::edges(graph.g); ep.first != ep.second; ++ep.first) {
        auto e = *ep.first;
        auto source = boost::source(e, graph.g);
        auto target = boost::target(e, graph.g);
        auto& edge_info = graph.g[e];

        QCEdge qc_edge = boost::add_edge(source, target, QCgraph).first;
        auto& qc_edge_info = QCgraph[qc_edge];
        qc_edge_info.c = edge_info.cost_derivative / 2;
        qc_edge_info.d = edge_info.cost;

        qc_edge_info.d -= edge_info.cost_derivative * od_set.link_flows[origin](source, target);
    }
}

template<typename cost_type>
double MQCF<cost_type>::compute_initial_delta() {
    double max_excess = QCgraph[origin].excess;
    int n = boost::num_vertices(QCgraph);
    int m = boost::num_edges(QCgraph);
    return max_excess / (2 * n + m);
}

template<typename cost_type>
void MQCF<cost_type>::adjust_flow(const double& delta) {
    auto edges = boost::edges(QCgraph);
    for (auto eit = edges.first; eit != edges.second; ++eit) {
        QCEdge e = *eit;
        QCVertex u = source(e, QCgraph);
        QCVertex v = target(e, QCgraph);
        double C_prime = 2 * QCgraph[e].c * (QCgraph[e].flow + delta) + QCgraph[e].d;

        if (C_prime < QCgraph[v].pi - QCgraph[u].pi) {
            QCgraph[e].flow += delta;
            QCgraph[u].excess -= delta;
            QCgraph[v].excess += delta;
        } else if (QCgraph[e].flow >= delta &&
                   (2 * QCgraph[e].c * (QCgraph[e].flow - delta) + QCgraph[e].d > QCgraph[v].pi - QCgraph[u].pi)) {
            QCgraph[e].flow -= delta;
            QCgraph[u].excess += delta;
            QCgraph[v].excess -= delta;
        }
        // Else keep flow unchanged
    }
}

template<typename cost_type>
void MQCF<cost_type>::augment_flow(const std::vector<std::pair<QCEdge, bool>>& pred, double delta, QCVertex& s, QCVertex& t) {
    QCVertex cur = t;

    while (cur != s) {
        QCEdge e = pred[cur].first;
        bool is_inverse = pred[cur].second;
        QCVertex prev = is_inverse ? target(e, QCgraph) : source(e, QCgraph);

        if (is_inverse) {
            // QCgraph[orig].reverse_flow += delta;
            QCgraph[e].flow -= delta;

            QCgraph[source(e, QCgraph)].excess += delta;
            QCgraph[target(e, QCgraph)].excess -= delta;
        } else {
            // Forward edge: increase flow
            QCgraph[e].flow += delta;

            QCgraph[source(e, QCgraph)].excess -= delta;
            QCgraph[target(e, QCgraph)].excess += delta;
        }

        cur = prev;
    }
}

template<typename cost_type>
void MQCF<cost_type>::update_potentials(const std::vector<double>& dist) {
    auto vertices = boost::vertices(QCgraph);
    for (auto vit = vertices.first; vit != vertices.second; ++vit) {
        QCgraph[*vit].pi = dist[*vit];  // Directly use shortest path distances as potentials
    }
}

template<typename cost_type>
void MQCF<cost_type>::basic_algorithm(int max_iter, double epsilon) {
    double Delta = compute_initial_delta();

    int n = boost::num_vertices(QCgraph);
    int m = boost::num_edges(QCgraph);

    int iteration = 0;

    while (Delta > epsilon / (2 * n + m + 1) && iteration < max_iter) {
        iteration += 1;
        // Main phase loop
        while (true) {
            // Find S(Δ) and T(Δ)
            std::vector<QCVertex> S, T;
            auto vertices = boost::vertices(QCgraph);
            for (auto vit = vertices.first; vit != vertices.second; ++vit) {
                if (QCgraph[*vit].excess >= Delta) S.push_back(*vit);
                if (QCgraph[*vit].excess <= -Delta) T.push_back(*vit);
            }

            if (S.empty() || T.empty()) break;

            // Find shortest path using Bellman-Ford (simple approach)
            std::vector<double> dist(n, INFINITY);
            std::vector<std::pair<QCEdge, bool>> pred(n);
            dist[S[0]] = 0;
            std::vector<QCVertex> node_pred(n);

            auto edges = boost::edges(QCgraph);
            for (size_t i = 0; i < n; ++i) {
                for (auto eit = edges.first; eit != edges.second; ++eit) {
                    QCEdge e = *eit;
                    auto& edge_info = QCgraph[e];
                    QCVertex u = source(e, QCgraph);
                    QCVertex v = target(e, QCgraph);
                    double cost = 2 * edge_info.c * (edge_info.flow + Delta) + edge_info.d;
                    double reverse_cost = 2 * edge_info.c * (edge_info.flow - Delta) + edge_info.d;
                    if (dist[u] + cost < dist[v] - 1e-10) {
                        dist[v] = dist[u] + cost;
                        pred[v] = {e, false};
                        node_pred[v] = u;
                    }
                    if (edge_info.flow >= Delta){
                        if (dist[v] - reverse_cost < dist[u] - 1e-10) {
                            dist[u] = dist[v] - reverse_cost;
                            pred[u] = {e, true};
                            node_pred[u] = v;
                        }
                    }
                }
            }

            update_potentials(dist);
            if (dist[T[0]] != INFINITY) {
                augment_flow(pred, Delta, S[0], T[0]);
            }
        }

        std::vector<double> dist(n, INFINITY);
        std::vector<std::pair<QCEdge, bool>> pred(n);
        dist[origin] = 0;
        std::vector<QCVertex> node_pred(n);
        auto edges = boost::edges(QCgraph);

        for (size_t i = 0; i < n; ++i) {
            for (auto eit = edges.first; eit != edges.second; ++eit) {
                QCEdge e = *eit;
                auto& edge_info = QCgraph[e];
                QCVertex u = source(e, QCgraph);
                QCVertex v = target(e, QCgraph);
                double cost = 2 * edge_info.c * (edge_info.flow + Delta) + edge_info.d;
                double reverse_cost = 2 * edge_info.c * (edge_info.flow - Delta) + edge_info.d;
                if (dist[u] + cost < dist[v]) {
                    dist[v] = dist[u] + cost;
                    pred[v] = {e, false};
                    node_pred[v] = u;
                }
                if (edge_info.flow >= Delta){
                    if (dist[v] - reverse_cost < dist[u]) {
                        dist[u] = dist[v] - reverse_cost;
                        pred[u] = {e, true};
                        node_pred[u] = v;
                    }
                }
            }
        }
        update_potentials(dist);  // Update node potentials

        adjust_flow(Delta / 2);
        Delta /= 2;
    }
    // print_edges();
}

template<typename cost_type> void MQCF<cost_type>::print_vertices() {
    std::cout << std::left << std::setw(10) << "ID" << std::left << std::setw(10) << "Excess" << std::endl;
    auto vertices = boost::vertices(QCgraph);
    for (auto it = vertices.first; it != vertices.second; ++it) {
        std::cout << std::left << std::setw(10) << *it << std::left << std::setw(10) << QCgraph[*it].excess << std::endl;
    }
}

template<typename cost_type> void MQCF<cost_type>::print_edges() {
    std::cout << std::left << std::setw(10) << "Source"
              << std::left << std::setw(10) << "Target"
              << std::left << std::setw(10) << "Flow" << std::endl;
    auto edges = boost::edges(QCgraph);
    for (auto it = edges.first; it != edges.second; ++it) {
        auto source = boost::source(*it, QCgraph);
        auto target = boost::target(*it, QCgraph);
        const auto& edge_info = QCgraph[*it];
        std::cout << std::left << std::setw(10) << source
                  << std::left << std::setw(10) << target
                  << std::left << std::setw(10) << edge_info.flow << std::endl;
    }
}

#endif //MQCF_H
