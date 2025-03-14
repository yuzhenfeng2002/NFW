//
// Created by Yuzhen Feng on 19/2/2025.
//

#ifndef IO_H
#define IO_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <regex>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "graph.h"
#include "demand.h"

template<typename cost_type>
void read_graph_from_file(const std::string& filename, Graph<cost_type>& graph) {
    // Check if the filename matches the pattern xxx_net.tntp
    if (std::regex pattern(".*_net\\.tntp"); !std::regex_match(filename, pattern)) {
        throw std::invalid_argument("Filename does not satisfy the pattern xxx_net.tntp");
    }

    // Open the file
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    int num_vertices = 0;
    int num_edges = 0;

    // Read the number of nodes
    while (std::getline(file, line)) {
        if (line.find("<NUMBER OF NODES>") != std::string::npos) {
            std::istringstream ss(line);
            std::string temp;
            ss >> temp >> temp >> temp >> num_vertices;
            break;
        }
    }
    graph.num_vertices = num_vertices;

    // Add nodes to the graph
    for (int i = 0; i < num_vertices; ++i) {
        boost::add_vertex(graph.g);
    }

    // Skip metadata lines
    while (std::getline(file, line)) {
        if (line.find("<END OF METADATA>") != std::string::npos) {
            break;
        }
    }

    // Read the graph data
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '~') {
            continue;
        }

        std::istringstream ss(line);
        int init_node, term_node;
        double capacity, length, fft, b, power, speed, toll;
        int link_type;

        ss >> init_node >> term_node >> capacity >> length >> fft >> b >> power >> speed >> toll >> link_type;

        int init_node_id = init_node - 1;
        int term_node_id = term_node - 1;

        // Add edge to the graph
        typename Graph<cost_type>::edge_type e;
        bool inserted;
        boost::tie(e, inserted) = boost::add_edge(
            boost::vertex(init_node_id, graph.g), boost::vertex(term_node_id, graph.g), graph.g);
        if (inserted) {
            graph.g[e].cost_fun.initialize(capacity, fft, b, power);
            graph.g[e].flow = 0.0;
            graph.g[e].update_flow(0.0);
            num_edges++;
        }
    }
    graph.num_edges = num_edges;
}

template<typename cost_type>
void read_demand_from_file(const std::string& filename,
    Graph<cost_type>& graph, OD_set<cost_type>& od_set) {
    // Check if the filename matches the pattern xxx_trips.tntp
    if (std::regex pattern(".*_trips\\.tntp"); !std::regex_match(filename, pattern)) {
        throw std::invalid_argument("Filename does not satisfy the pattern xxx_trips.tntp");
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;

    // Skip metadata lines
    while (std::getline(file, line)) {
        if (line.find("<END OF METADATA>") != std::string::npos) {
            break;
        }
    }

    // Read the demand data
    int origin_id = -1;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '~') {
            continue;
        }

        if (line.find("Origin") != std::string::npos) {
            std::istringstream iss(line);
            std::string temp;
            iss >> temp >> origin_id;
            origin_id -= 1; // Adjust for zero-based indexing
        } else {
            std::istringstream iss(line);
            std::string app;
            while (getline(iss, app, ';')) {
                boost::trim(app);
                if (app.empty()) {
                    continue;
                }

                std::vector<std::string> result;
                boost::split(result, app, boost::is_any_of(":"), boost::token_compress_on);

                boost::trim(result[0]);
                boost::trim(result[1]);
                int destination_id = boost::lexical_cast<int>(result[0]) - 1;
                auto flow = boost::lexical_cast<double>(result[1]);

                if (origin_id != destination_id && flow > 0.0) {
                    typename  Graph<cost_type>::vertex_type origin = boost::vertex(origin_id, graph.g);
                    typename  Graph<cost_type>::vertex_type destination = boost::vertex(destination_id, graph.g);
                    od_set.add_od_pair(origin, destination, flow);
                }
            }
        }
    }

    int last_origin_id = od_set.od_pairs[0].origin;
    int i = 0;
    od_set.ods_from_origin.push_back({});
    for (auto& od : od_set.od_pairs) {
        if (od.origin != last_origin_id) {
            last_origin_id = od.origin;
            od_set.ods_from_origin.push_back({});
            i++;
        }
        od_set.ods_from_origin[i].push_back(&od);
    }
}

#endif //IO_H
