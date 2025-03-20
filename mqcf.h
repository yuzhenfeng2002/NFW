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
const double INF = std::numeric_limits<double>::max();

struct NodeProps {
    double demand{0};      // Node demand (b_i)
    double excess{0};      // Current excess (ρ_f(i) - b_i)
    double pi{0};          // Dual variable (π_i)
    bool is_S{false};
    bool is_T{false};
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
    int m, n;
    std::pair<boost::graph_traits<QCGraph>::edge_iterator, boost::graph_traits<QCGraph>::edge_iterator> qc_edges;

    MQCF() : origin(0), m(0), n(0) {}
    MQCF(const Graph<cost_type>& graph, const OD_set<cost_type>& od_set, int origin, int idx);
    void update_graph(const Graph<cost_type>& graph, const OD_set<cost_type>& od_set, int origin, int idx);

    void basic_algorithm(int max_iter=100, double epsilon=1e-6);
private:
    double compute_initial_delta();
    void dijkstra_for_MQCF(const QCVertex& s, QCVertex& t, std::vector<double>& dist, std::vector<std::pair<QCEdge, bool>>& pred, double Delta);
    void revised_dijkstra_for_MQCF(const std::vector<QCVertex>& S, QCVertex& s, QCVertex& t, std::vector<double>& dist, std::vector<std::pair<QCEdge, bool>>& pred, double Delta);
    void adjust_flow(const double& delta);
    void augment_flow(const std::vector<std::pair<QCEdge, bool>>& pred, double delta, QCVertex& s, QCVertex& t);
    void update_potentials(const std::vector<double>& dist);
    bool check_potentials(const double& delta);
    void print_edges();
    void print_vertices();
};

template<typename cost_type>
MQCF<cost_type>::MQCF(const Graph<cost_type>& graph, const OD_set<cost_type>& od_set, int origin, int idx) {
    this->origin = origin;
    n = boost::num_vertices(graph.g);
    m = boost::num_edges(graph.g);

    for (auto vp = boost::vertices(graph.g); vp.first != vp.second; ++vp.first) {
        auto v = *vp.first;
        boost::add_vertex(QCgraph);
    }

    for (auto od : od_set.ods_from_origin[idx]) {
        auto destination = od->destination;
        auto flow = od->flow;
        QCgraph[destination].demand += flow;
        QCgraph[origin].demand -= flow;
    }

    for (auto ep = boost::edges(graph.g); ep.first != ep.second; ++ep.first) {
        auto e = *ep.first;
        auto source = boost::source(e, graph.g);
        auto target = boost::target(e, graph.g);
        boost::add_edge(source, target, QCgraph);
    }

    qc_edges = boost::edges(QCgraph);
}

template<typename cost_type>
void MQCF<cost_type>::update_graph(const Graph<cost_type>& graph, const OD_set<cost_type>& od_set, int origin, int idx) {
    QCgraph[origin].excess = - QCgraph[origin].demand;
    for (auto od : od_set.ods_from_origin[idx]) {
        auto destination = od->destination;
        QCgraph[destination].excess = - QCgraph[destination].demand;
    }
    qc_edges = boost::edges(QCgraph);
    for (auto ep = qc_edges; ep.first != ep.second; ++ep.first) {
        auto e = *ep.first;
        auto source = boost::source(e, graph.g);
        auto target = boost::target(e, graph.g);
        auto& qc_edge_info = QCgraph[e];

        auto edge = boost::edge(source, target, graph.g).first;
        auto& edge_info = graph.g[edge];
        qc_edge_info.c = edge_info.cost_derivative / 2;
        qc_edge_info.d = edge_info.cost;
        qc_edge_info.d -= edge_info.cost_derivative * edge_info.orgin_flows[idx];
        qc_edge_info.flow = 0;
    }
}

template<typename cost_type>
void MQCF<cost_type>::adjust_flow(const double& delta) {
    for (auto eit = qc_edges.first; eit != qc_edges.second; ++eit) {
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
        } else {
            // Forward edge: increase flow
            QCgraph[e].flow += delta;
        }
        cur = prev;
    }
    QCgraph[s].excess -= delta;
    if (QCgraph[s].excess < delta) {
        QCgraph[s].is_S = false;
    }
    QCgraph[t].excess += delta;
    if (QCgraph[t].excess > -delta) {
        QCgraph[t].is_T = false;
    }
}

template<typename cost_type>
void MQCF<cost_type>::update_potentials(const std::vector<double>& dist) {
    auto vertices = boost::vertices(QCgraph);
    for (auto vit = vertices.first; vit != vertices.second; ++vit) {
        QCgraph[*vit].pi = dist[*vit];  // Directly use shortest path distances as potentials
    }
}

inline std::pair<bool, std::vector<QCEdge>> bellman_ford_for_MQCF(const QCGraph& QCgraph, const QCVertex& origin, std::vector<double>& dist, std::vector<std::pair<QCEdge, bool>>& pred, double Delta) {
    int n = boost::num_vertices(QCgraph);

    dist.assign(n, INF);
    dist[origin] = 0;

    auto edges = boost::edges(QCgraph);
    for (size_t i = 0; i < n - 1; ++i) {
        bool updated = false;
        for (auto eit = edges.first; eit != edges.second; ++eit) {
            QCEdge e = *eit;
            auto& edge_info = QCgraph[e];
            QCVertex u = source(e, QCgraph);
            QCVertex v = target(e, QCgraph);
            double cost = 2 * edge_info.c * (edge_info.flow + Delta) + edge_info.d;
            if (dist[u] + cost < dist[v] - 1e-10) {
                dist[v] = dist[u] + cost;
                pred[v] = {e, false};
                updated = true;
            }
        }
        if (!updated) return {false, {}};
    }

    for (auto eit = edges.first; eit != edges.second; ++eit) {
        QCEdge e = *eit;
        auto& edge_info = QCgraph[e];
        QCVertex u = source(e, QCgraph);
        QCVertex v = target(e, QCgraph);
        double cost = 2 * edge_info.c * (edge_info.flow + Delta) + edge_info.d;
        if (dist[u] != std::numeric_limits<double>::max() && dist[u] + cost < dist[v]) {
            std::vector<QCEdge> cycle;
            QCVertex x = v;
            for (int i = 0; i < n; ++i) {
                x = pred[x].first.m_source;
            }
            for (QCVertex cur = x; ; cur = pred[cur].first.m_source) {
                QCEdge pred_edge = pred[cur].first;
                if (cur == x && cycle.size() > 1) break;
                cycle.push_back(pred_edge);
            }
            std::reverse(cycle.begin(), cycle.end());
            return std::make_pair(true, cycle);
        }
    }
    return {false, {}};
}

template<typename cost_type>
void MQCF<cost_type>::dijkstra_for_MQCF(const QCVertex& s, QCVertex& t, std::vector<double>& dist, std::vector<std::pair<QCEdge, bool>>& pred, double Delta) {
    typedef std::pair<double, QCVertex> VertexCostPair;
    std::vector<bool> visited(boost::num_vertices(QCgraph), false);

    std::vector<double> modified_dist(boost::num_vertices(QCgraph), INF);
    modified_dist[s] = 0;
    dist[s] = 0;
    std::priority_queue<VertexCostPair, std::vector<VertexCostPair>, std::greater<VertexCostPair>> pq;
    pq.push({0, s});

    while (!pq.empty()) {
        QCVertex u = pq.top().second;
        pq.pop();
        if (visited[u]) continue;

        /* Bug should be fixed here */
        // if (QCgraph[u].is_T) {
        //     t = u;
        //     auto vertices = boost::vertices(QCgraph);
        //     double alpha = dist[u] - QCgraph[u].pi;
        //     while (!pq.empty()) {
        //         QCVertex v = pq.top().second;
        //         pq.pop();
        //         if (dist[v] < INF) {
        //             alpha = std::min(alpha, dist[v] - QCgraph[v].pi);
        //         }
        //     }
        //     for (auto vit = vertices.first; vit != vertices.second; ++vit) {
        //         if (!visited[*vit]) {
        //             QCgraph[*vit].pi += alpha;
        //         }
        //     }
        //     break;
        // }

        visited[u] = true;
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
                pq.push({modified_dist[u] + modified_cost, v});
            }
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
                pq.push({modified_dist[u] + modified_cost, v});
            }
        }

        QCgraph[u].pi = dist[u];
    }
    QCgraph[s].pi = dist[s];
    // if (!check_potentials(Delta)) throw std::runtime_error("Potentials are not consistent!");
}

template<typename cost_type>
void MQCF<cost_type>::revised_dijkstra_for_MQCF(const std::vector<QCVertex> &S, QCVertex &s, QCVertex &t,
    std::vector<double> &dist, std::vector<std::pair<QCEdge, bool>> &pred, double Delta) {
    typedef std::pair<double, std::pair<QCEdge, bool>> EdgeCostPair;
    std::vector<bool> visited(boost::num_vertices(QCgraph), false);
    std::list<EdgeCostPair> pq;
    for (auto u : S) {
        if (QCgraph[u].is_S) {
            visited[u] = true;
        }
    }
    for (auto u : S) {
        if (!QCgraph[u].is_S) continue;
        auto out_edges = boost::out_edges(u, QCgraph);
        auto in_edges = boost::in_edges(u, QCgraph);
        for (auto eit = out_edges.first; eit != out_edges.second; ++eit) {
            QCEdge e = *eit;
            QCVertex v = boost::target(e, QCgraph);
            if (visited[v]) continue;
            const auto& edge_info = QCgraph[e];
            double cost = 2 * edge_info.c * (edge_info.flow + Delta) + edge_info.d;
            double modified_cost = cost - QCgraph[v].pi + QCgraph[u].pi;
            pq.push_back({modified_cost, {e, false}});
        }
        for (auto eit = in_edges.first; eit != in_edges.second; ++eit) {
            QCEdge e = *eit;
            QCVertex v = boost::source(e, QCgraph);
            if (visited[v]) continue;
            const auto& edge_info = QCgraph[e];
            if (edge_info.flow < Delta) continue;
            double cost = - 2 * edge_info.c * (edge_info.flow - Delta) - edge_info.d;
            double modified_cost = cost - QCgraph[v].pi + QCgraph[u].pi;
            pq.push_back({modified_cost, {e, true}});
        }
    }
    while (!pq.empty()) {
        /* For debugging */
        // auto vertices = boost::vertices(QCgraph);
        // auto edges = boost::edges(QCgraph);
        // std::cout << "Vertex\tVisited\tExcess\tPi" << std::endl;
        // for (auto vit = vertices.first; vit != vertices.second; ++vit) {
        //     std::cout << *vit << '\t' << visited[*vit] << '\t' << QCgraph[*vit].excess << '\t' << QCgraph[*vit].pi << std::endl;
        // }
        // std::cout << "Edge\tFlow\tCost" << std::endl;
        // for (auto eit = edges.first; eit != edges.second; ++eit) {
        //     QCEdge e = *eit;
        //     auto source = boost::source(e, QCgraph);
        //     auto target = boost::target(e, QCgraph);
        //     auto& edge_info = QCgraph[e];
        //     double cost = 2 * edge_info.c * (edge_info.flow + Delta) + edge_info.d;
        //     std::cout << source << '\t' << target << '\t' << edge_info.flow << '\t' << cost << std::endl;
        //     if (edge_info.flow > Delta) {
        //         cost = - 2 * edge_info.c * (edge_info.flow - Delta) - edge_info.d;
        //         std::cout << target << '\t' << source << '\t' << edge_info.flow << '\t' << cost << std::endl;
        //     }
        // }


        double min_cost = std::numeric_limits<double>::max();
        auto min_it = pq.end();
        for (auto it = pq.begin(); it != pq.end(); ++it) {
            if (it->first < min_cost) {
                min_cost = it->first;
                min_it = it;
            }
        }
        auto alpha = min_it->first;
        auto e = min_it->second.first;
        auto is_reverse = min_it->second.second;
        pq.erase(min_it);
        QCVertex u, v;
        if (is_reverse) {
            v = boost::source(e, QCgraph);
            u = boost::target(e, QCgraph);
        } else {
            u = boost::source(e, QCgraph);
            v = boost::target(e, QCgraph);
        }
        if (visited[v]) continue;

        if (alpha > 0) {
            for (auto& edge_cost_pair : pq) {
                edge_cost_pair.first -= alpha;
            }
            auto vertices = boost::vertices(QCgraph);
            for (auto vit = vertices.first; vit != vertices.second; ++vit) {
                if (!visited[*vit]) {
                    QCgraph[*vit].pi += alpha;
                }
            }
        }

        pred[v] = {e, is_reverse};

        if (QCgraph[v].is_T) {
            t = v;
            s = t;
            while (!QCgraph[s].is_S) {
                s = pred[s].second? boost::target(pred[s].first, QCgraph) : boost::source(pred[s].first, QCgraph);
            }
            break;
        }

        u = v;
        auto out_edges = boost::out_edges(u, QCgraph);
        auto in_edges = boost::in_edges(u, QCgraph);
        for (auto eit = out_edges.first; eit != out_edges.second; ++eit) {
            QCEdge e = *eit;
            QCVertex v = boost::target(e, QCgraph);
            if (visited[v]) continue;
            const auto& edge_info = QCgraph[e];
            double cost = 2 * edge_info.c * (edge_info.flow + Delta) + edge_info.d;
            double modified_cost = cost - QCgraph[v].pi + QCgraph[u].pi;
            pq.push_back({modified_cost, {e, false}});
        }
        for (auto eit = in_edges.first; eit != in_edges.second; ++eit) {
            QCEdge e = *eit;
            QCVertex v = boost::source(e, QCgraph);
            if (visited[v]) continue;
            const auto& edge_info = QCgraph[e];
            if (edge_info.flow < Delta) continue;
            double cost = - 2 * edge_info.c * (edge_info.flow - Delta) - edge_info.d;
            double modified_cost = cost - QCgraph[v].pi + QCgraph[u].pi;
            pq.push_back({modified_cost, {e, true}});
        }
        visited[v] = true;
    }
}

template<typename cost_type>
double MQCF<cost_type>::compute_initial_delta() {
    double max_excess = QCgraph[origin].excess;
    double Delta = max_excess / (2 * n + m);
    std::vector<double> dist(n, INF);
    std::vector<std::pair<QCEdge, bool>> pred(n);
    do {
        auto [has_cycle, cycle] = bellman_ford_for_MQCF(QCgraph, origin, dist, pred, Delta);
        if (!has_cycle) break;
        else {
            double two_c = 0, d = 0;
            for (auto& edge : cycle) {
                auto& edge_info = QCgraph[edge];
                two_c += edge_info.c * 2;
                d += edge_info.d;
            }
            Delta = std::max(Delta, -d / two_c);
        }
    } while (true);
    update_potentials(dist);
    return Delta;
}

template<typename cost_type>
void MQCF<cost_type>::basic_algorithm(int max_iter, double epsilon) {
    double Delta = compute_initial_delta();
    int iteration = 0;
    std::vector<double> dist(n, INF);
    std::vector<std::pair<QCEdge, bool>> pred(n);
    auto vertices = boost::vertices(QCgraph);
    while (Delta > epsilon / (2 * n + m + 1) && iteration < max_iter) {
        iteration += 1;
        // Main phase loop
        std::vector<QCVertex> S, T;
        for (auto vit = vertices.first; vit != vertices.second; ++vit) {
            auto& vertex_info = QCgraph[*vit];
            if (vertex_info.excess >= Delta) {
                vertex_info.is_S = true;
                S.push_back(*vit);
            } else vertex_info.is_S = false;
            if (vertex_info.excess <= -Delta) {
                vertex_info.is_T = true;
                T.push_back(*vit);
            } else vertex_info.is_T = false;
        }
        while (true) {
            // Find S(Δ) and T(Δ)
            if (S.empty() || T.empty()) break;
            while (!S.empty() && !QCgraph[S[0]].is_S) {
                S.erase(S.begin());
            }
            if (S.empty()) break;
            while (!T.empty() && !QCgraph[T[0]].is_T) {
                T.erase(T.begin());
            }
            if (T.empty()) break;
            // std::vector<int> pred_node(n);
            // bellman_ford_for_MQCF(QCgraph, S[0], dist, pred, Delta);
            // for (size_t i = 1; i < n; ++i) {
            //     pred_node[i] = pred[i].first.m_source == i ? pred[i].first.m_target : pred[i].first.m_source;
            // }
            // std::vector<int> tmp_pred_node(n);
            // std::vector<double> tmp(n, INFINITY);
            QCVertex s = S[0];
            QCVertex t = T[0];
            // dijkstra_for_MQCF(s, t, dist, pred, Delta);
            revised_dijkstra_for_MQCF(S, s, t, dist, pred, Delta);
            // for (size_t i = 1; i < n; ++i) {
            //     tmp_pred_node[i] = pred[i].first.m_source == i ? pred[i].first.m_target : pred[i].first.m_source;
            // }
            // if (!check_potentials(Delta)) throw std::runtime_error("Potentials are not consistent!");
            augment_flow(pred, Delta, s, t);

        }
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

template<typename cost_type> bool MQCF<cost_type>::check_potentials(const double& delta) {
    for (auto eit = qc_edges.first; eit != qc_edges.second; ++eit) {
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
