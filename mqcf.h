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

const double EPS_NUM = 1e-2;

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

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, NodeProps, EdgeProps> QCGraph;
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
    bool check_potentials(const double& delta);
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

        qc_edge_info.d -= edge_info.cost_derivative * edge_info.orgin_flows[origin];
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

inline void bellman_ford_for_MQCF(const QCGraph& QCgraph, const QCVertex& origin, std::vector<double>& dist, std::vector<std::pair<QCEdge, bool>>& pred, double Delta) {
    int n = boost::num_vertices(QCgraph);

    const double inf = std::numeric_limits<double>::max();
    dist.assign(num_vertices(QCgraph), inf);
    dist[origin] = 0;

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
            }
            if (edge_info.flow >= Delta){
                if (dist[v] - reverse_cost < dist[u] - 1e-10) {
                    dist[u] = dist[v] - reverse_cost;
                    pred[u] = {e, true};
                }
            }
        }
    }
}

inline void dijkstra_for_MQCF(const QCGraph& QCgraph, const QCVertex& origin, std::vector<double>& dist, std::vector<std::pair<QCEdge, bool>>& pred, double Delta) {
    typedef std::pair<double, QCVertex> VertexCostPair;
    std::vector<bool> visited(boost::num_vertices(QCgraph), false);
    const double inf = std::numeric_limits<double>::max();

    std::vector<double> modified_dist(boost::num_vertices(QCgraph), inf);
    modified_dist.assign(num_vertices(QCgraph), inf);
    modified_dist[origin] = 0;
    dist[origin] = 0;
    visited[origin] = true;
    std::priority_queue<VertexCostPair, std::vector<VertexCostPair>, std::greater<VertexCostPair>> pq;
    pq.push({0, origin});

    while (!pq.empty()) {
        QCVertex u = pq.top().second;
        visited[u] = true;
        pq.pop();

        auto out_edges = boost::out_edges(u, QCgraph);
        for (auto eit = out_edges.first; eit != out_edges.second; ++eit) {
            QCEdge e = *eit;
            QCVertex v = boost::target(e, QCgraph);
            if (visited[v]) continue;
            const auto& edge_info = QCgraph[e];
            double cost = 2 * edge_info.c * (edge_info.flow + Delta) + edge_info.d;
            double modified_cost = cost - QCgraph[v].pi + QCgraph[u].pi;

            if (modified_dist[u] + modified_cost < modified_dist[v]) {
                modified_dist[v] = modified_dist[u] + modified_cost;
                dist[v] = dist[u] + cost;
                pred[v] = {e, false};
            }
            pq.push({modified_dist[u] + modified_cost, v});
        }

        auto in_edges = boost::in_edges(u, QCgraph);
        for (auto eit = in_edges.first; eit != in_edges.second; ++eit) {
            QCEdge e = *eit;
            QCVertex v = boost::source(e, QCgraph);
            if (visited[v]) continue;
            const auto& edge_info = QCgraph[e];
            if (edge_info.flow < Delta) continue;
            double cost = - 2 * edge_info.c * (edge_info.flow - Delta) - edge_info.d;
            double modified_cost = cost - QCgraph[v].pi + QCgraph[u].pi;

            if (modified_dist[u] + modified_cost < modified_dist[v]) {
                modified_dist[v] = modified_dist[u] + modified_cost;
                dist[v] = dist[u] + cost;
                pred[v] = {e, true};
            }
            pq.push({modified_dist[u] + modified_cost, v});
        }
    }
}

template<typename cost_type>
void MQCF<cost_type>::basic_algorithm(int max_iter, double epsilon) {
    double Delta = compute_initial_delta() * 32;

    int n = boost::num_vertices(QCgraph);
    int m = boost::num_edges(QCgraph);

    int iteration = 0;

    std::vector<double> dist(n, INFINITY);
    std::vector<std::pair<QCEdge, bool>> pred(n);

    bellman_ford_for_MQCF(QCgraph, origin, dist, pred, Delta);
    update_potentials(dist);
    if (!check_potentials(Delta)) {
        throw std::runtime_error("Potential check failed");
    }

    while (Delta > epsilon / (2 * n + m + 1) && iteration < max_iter) {
        iteration += 1;
        // Main phase loop
        int k = 0;
        while (true) {
            // Find S(Δ) and T(Δ)
            std::vector<QCVertex> S, T;
            auto vertices = boost::vertices(QCgraph);
            for (auto vit = vertices.first; vit != vertices.second; ++vit) {
                if (QCgraph[*vit].excess >= Delta) S.push_back(*vit);
                if (QCgraph[*vit].excess <= -Delta) T.push_back(*vit);
            }

            if (S.empty() || T.empty()) break;

            // std::vector<int> pred_node(n);
            // bellman_ford_for_MQCF(QCgraph, S[0], dist, pred, Delta);
            // for (size_t i = 1; i < n; ++i) {
            //     pred_node[i] = pred[i].first.m_source == i ? pred[i].first.m_target : pred[i].first.m_source;
            // }
            // std::vector<int> tmp_pred_node(n);
            // std::vector<double> tmp(n, INFINITY);
            dijkstra_for_MQCF(QCgraph, S[0], dist, pred, Delta);
            // for (size_t i = 1; i < n; ++i) {
            //     tmp_pred_node[i] = pred[i].first.m_source == i ? pred[i].first.m_target : pred[i].first.m_source;
            // }
            update_potentials(dist);

            if (dist[T[0]] != INFINITY) {
                augment_flow(pred, Delta, S[0], T[0]);
            }
            if (!check_potentials(Delta)) {
                throw std::runtime_error("Potential check failed");
            }
            k++;
        }

        dijkstra_for_MQCF(QCgraph, origin, dist, pred, Delta);
        update_potentials(dist);

        adjust_flow(Delta / 2);
        Delta /= 2;
        if (!check_potentials(Delta)) {
            throw std::runtime_error("Potential check failed");
        }
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

template<typename cost_type> bool MQCF<cost_type>::check_potentials(const double& delta) {
    auto edges = boost::edges(QCgraph);
    for (auto eit = edges.first; eit != edges.second; ++eit) {
        QCEdge e = *eit;
        QCVertex u = source(e, QCgraph);
        QCVertex v = target(e, QCgraph);
        auto& edge_info = QCgraph[e];
        double cost = 2 * edge_info.c * (edge_info.flow + delta) + edge_info.d;
        double reverse_cost = 2 * edge_info.c * (edge_info.flow - delta) + edge_info.d;
        if (cost < QCgraph[v].pi - QCgraph[u].pi - EPS_NUM) {
            std::cout << '<' << u << ", " << v << "> cost: " << cost << " pi_u: " << QCgraph[u].pi << " pi_v: " << QCgraph[v].pi << std::endl;
            return false;
        }

        if (edge_info.flow >= delta && reverse_cost > QCgraph[v].pi - QCgraph[u].pi + EPS_NUM) {
            std::cout << '<' << v << ", " << u << "> cost: " << -reverse_cost << " pi_v: " << QCgraph[v].pi << " pi_u: " << QCgraph[u].pi << std::endl;
            return false;
        }
    }
    return true;
}

#endif //MQCF_H
