//
// Created by 封钰震 on 9/3/25.
//

#ifndef MQCF_H
#define MQCF_H

#include <iostream>
#include <vector>
#include <queue>
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
    bool is_reverse{false};    // Marks reverse edges
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, NodeProps, EdgeProps> QCGraph;
typedef typename boost::graph_traits<QCGraph>::vertex_descriptor QCVertex;
typedef typename boost::graph_traits<QCGraph>::edge_descriptor QCEdge;

template<typename cost_type>
class MQCF {
public:
    QCGraph QCgraph;
    int origin;
    std::vector<std::pair<int, double>> destinations;

    MQCF(const Graph<cost_type>& graph, const OD_set<cost_type>& od_set, int& origin);

    void basic_algorithm(int max_iter=20, double epsilon=1e-3);
private:
    double compute_initial_delta();
    void adjust_flow(const double& delta);
    void build_residual_graph(double Delta, std::vector<QCEdge>& residual_edges);
    void augment_flow(const std::vector<QCVertex>& pred, double delta, QCVertex& s, QCVertex& t);
    void update_potentials(const std::vector<double>& dist);
};

template<typename cost_type>
MQCF<cost_type>::MQCF(const Graph<cost_type>& graph, const OD_set<cost_type>& od_set, int& origin) {
    this->origin = origin;

    for (auto vp = boost::vertices(graph.g); vp.first != vp.second; ++vp.first) {
        auto v = *vp.first;
        boost::add_vertex(QCgraph);
    }
    for (auto ep = boost::edges(graph.g); ep.first != ep.second; ++ep.first) {
        auto e = *ep.first;
        auto source = boost::source(e, graph.g);
        auto target = boost::target(e, graph.g);
        auto& edge_info = graph.g[e];

        QCEdge qc_edge = boost::add_edge(source, target, QCgraph).first;
        QCgraph[qc_edge].c = edge_info.cost_derivative / 2; // Assuming cost is the quadratic coefficient
        QCgraph[qc_edge].d = edge_info.cost; // Set linear coefficient to 0 or appropriate value
        QCgraph[qc_edge].is_reverse = false;
    }

    for (auto od : od_set.ods_from_origin[origin]) {
        auto destination = od->destination;
        auto flow = od->flow;
        destinations.push_back(std::make_pair(destination, flow));
        QCgraph[destination].demand += flow;
        QCgraph[destination].excess -= flow;
        QCgraph[origin].excess += flow;
    }

    auto edges = boost::edges(graph.g);
    for (auto it = edges.first; it != edges.second; ++it) {
        auto edge = *it;
        auto& edge_info = graph.g[edge];
        auto source = boost::source(edge, graph.g);
        auto target = boost::target(edge, graph.g);
        auto qc_edge = boost::edge(source, target, QCgraph).first;
        auto& qc_edge_info = QCgraph[qc_edge];
        for (auto od : od_set.ods_from_origin[origin]) {
            qc_edge_info.d -= edge_info.cost_derivative * od->link_flows(source, target);
        }
    }
}

template<typename cost_type>
double MQCF<cost_type>::compute_initial_delta() {
    double max_excess = 0;
    for (auto& dest : destinations) {
        max_excess += dest.second;
    }
    int n = boost::num_vertices(QCgraph);
    int m = boost::num_edges(QCgraph);
    return max_excess / (2 * n + m); // Simplified heuristic
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
        } else if (QCgraph[e].flow >= delta &&
                   (2 * QCgraph[e].c * (QCgraph[e].flow - delta) + QCgraph[e].d > QCgraph[v].pi - QCgraph[u].pi)) {
            QCgraph[e].flow -= delta;
        }
        // Else keep flow unchanged
    }
}

// Complete flow augmentation implementation
template<typename cost_type>
void MQCF<cost_type>::augment_flow(const std::vector<QCVertex>& pred, double delta, QCVertex& s, QCVertex& t) {
    QCVertex cur = t;

    // Walk the path from t to s
    while (cur != s) {
        QCVertex prev = pred[cur];

        // Find edge between prev->cur
        auto er = boost::edge(prev, cur, QCgraph);
        if (!er.second) {
            throw std::runtime_error("Invalid predecessor path");
        }
        QCEdge e = er.first;

        if (QCgraph[e].is_reverse) {
            // Reverse edge: decrease original edge's flow
            QCEdge orig = boost::edge(cur, prev, QCgraph).first;
            QCgraph[orig].flow -= delta;

            // Update excess for original nodes
            QCgraph[source(orig, QCgraph)].excess += delta;
            QCgraph[target(orig, QCgraph)].excess -= delta;
        } else {
            // Forward edge: increase flow
            QCgraph[e].flow += delta;

            // Update excess for endpoints
            QCgraph[source(e, QCgraph)].excess -= delta;
            QCgraph[target(e, QCgraph)].excess += delta;
        }

        cur = prev;
    }
}

// Updated residual network construction
template<typename cost_type>
void MQCF<cost_type>::build_residual_graph(double Delta, std::vector<QCEdge>& residual_edges) {
    residual_edges.clear();
    auto edges = boost::edges(QCgraph);

    for (auto eit = edges.first; eit != edges.second; ++eit) {
        QCEdge e = *eit;

        // Add forward edge if in original graph
        if (!QCgraph[e].is_reverse) {
            residual_edges.push_back(e);

            // Add reverse edge if residual capacity >= Δ
            if (QCgraph[e].flow >= Delta) {
                QCEdge reverse_e = boost::add_edge(target(e, QCgraph), source(e, QCgraph), QCgraph).first;
                QCgraph[reverse_e].is_reverse = true;
                QCgraph[reverse_e].c = QCgraph[e].c;
                QCgraph[reverse_e].d = -QCgraph[e].d;

                residual_edges.push_back(reverse_e);
            }
        }
    }
}

// Enhanced potential update
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
        std::vector<QCEdge> residual_edges;
        while (true) {
            build_residual_graph(Delta, residual_edges);
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
            std::vector<QCVertex> pred(n);
            dist[S[0]] = 0;

            // Relax edges (Bellman-Ford)
            for (size_t i = 0; i < n; ++i) {
                for (auto& e : residual_edges) {
                    QCVertex u = source(e, QCgraph);
                    QCVertex v = target(e, QCgraph);
                    double cost = 2 * QCgraph[e].c * (QCgraph[e].flow + Delta) + QCgraph[e].d;
                    if (dist[u] + cost <= dist[v]) {
                        dist[v] = dist[u] + cost;
                        pred[v] = u;
                    }
                }
            }

            update_potentials(dist);  // Update node potentials

            // Augment flow along shortest path
            for (QCVertex t : T) {
                if (dist[t] != INFINITY) {
                    augment_flow(pred, Delta, S[0], t);
                }
            }
        }

        adjust_flow(Delta/2);
        Delta /= 2;
    }
}

#endif //MQCF_H
