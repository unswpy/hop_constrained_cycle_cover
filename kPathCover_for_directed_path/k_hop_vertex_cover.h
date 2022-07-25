#include "graph.h"
//in all the k hop vertex cover problem, we only use adjacency_list and adjacency_list_reverse, we do not store adjacency_list_double
bool delete_adjacency_list_with_specific_edge(vector<NODE_TYPE>* adjacency_list, NODE_TYPE node1, NODE_TYPE node2);
set<NODE_TYPE> k_hop_vertex_cover(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num);
bool find_k_path_from_node1(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result);
unordered_set<NODE_TYPE> k_hop_vertex_cover_heuristic(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num);
NODE_TYPE find_result_node_in_single_k_path(vector<NODE_TYPE> temp_k_path, int* find_times);
void test_k_hop_cover(int k, char* dataset, char* algorithm);
bool node_neccessary(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result);
unordered_set<NODE_TYPE> k_hop_vertex_cover_node_necessary(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num);
bool find_k_path_from_node1_enum_paths(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, set<vector<NODE_TYPE>>& result);
bool node_neccessary_enurate_all_paths(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result);
bool judge_cover_correctness(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num, unordered_set<NODE_TYPE>& test_cover);
bool find_k_path_from_node1_enum_paths_node_orders(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, set<vector<NODE_TYPE>>& result, unordered_set<NODE_TYPE>& check_nodes, vector<NODE_TYPE>& node_orders);
bool find_k_path_from_node1_node_orders(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result, unordered_set<NODE_TYPE>& check_nodes, vector<NODE_TYPE>& node_orders);
void find_dfs_order(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result, vector<bool>& check_nodes, vector<NODE_TYPE>& node_orders);
unordered_set<NODE_TYPE> k_hop_vertex_cover_node_necessary_dfs_order(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num);
void find_dfs_order_pre(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result, vector<bool>& check_nodes, vector<NODE_TYPE>& node_orders);
unordered_set<NODE_TYPE> k_hop_vertex_cover_node_necessary_dfs_pre_order(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num);
void find_bfs_order(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result, vector<bool>& check_nodes, vector<NODE_TYPE>& node_orders);
unordered_set<NODE_TYPE> k_hop_vertex_cover_node_necessary_bfs_order(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num);