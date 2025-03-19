//
// Created by Yuzhen Feng on 18/2/2025.
//

#include <chrono>
#include <iostream>

#include "ue_gp.h"
#include "ue_fw.h"
#include "ue_nfw.h"
#include "params.h"
#include "graph.h"
#include "demand.h"
#include "io.h"
#include "cost.h"

int main()
{
    Graph<bpr> graph;
    read_graph_from_file("data/" + NET_NAME + "_net.tntp", graph);
    OD_set<bpr> od_set;
    read_demand_from_file("data/" + NET_NAME + "_trips.tntp", graph, od_set);

    graph.print_vertices();
    graph.print_edges();

    auto initialization_start = std::chrono::high_resolution_clock::now();
    // UE_GP<bpr> ue(graph, od_set);
    // UE_FW<bpr> ue(graph, od_set);
    UE_NFW<bpr> ue(graph, od_set);
    ue.initialization();
    auto initialization_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> initialization_elapsed = initialization_end - initialization_start;

    auto equilibration_start = std::chrono::high_resolution_clock::now();
    // ue.gradient_projection(MAX_ITER, EPS_EQ, STEP_SIZE);
    // ue.frank_wolfe(MAX_ITER * 100, EPS_EQ);
    ue.newton_frank_wolfe(MAX_ITER, EPS_EQ);
    auto equilibration_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> equilibration_elapsed = equilibration_end - equilibration_start;

    std::cout << "Initialization time: " << initialization_elapsed.count() << "s" << std::endl;
    std::cout << "Equilibration time: " << equilibration_elapsed.count() << "s" << std::endl;

    graph.print_edges();
    // ue.print_average_cost();
    return 0;
}
