#include "graph.h"
#include "BFS_prune.h"
int SHORT_CUT_DIS = 4;


int bfs_return_spd(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  k, NODE_TYPE query_node1, NODE_TYPE query_node2)
{//store the final induced subgraph into reverse_adjacency_in_subgraph_left
	int result = -1;
	DISTANCE_TYPE cur_distance = 0;
	set<NODE_TYPE> left_proprogate;
	set<NODE_TYPE> right_proprogate;
	//cur_proprogate.insert(query_node1);
	set<NODE_TYPE> node1_set;
	set<NODE_TYPE> node2_set;
	//set<NODE_TYPE> visited_nodes;

	set<NODE_TYPE> left_visited_nodes;
	set<NODE_TYPE> right_visited_nodes;

	//dag_min_induced_subgraph_reverse.insert(std::make_pair(query_node2, node2_set));// query node1's parent is empty
	left_visited_nodes.insert(query_node1);
	right_visited_nodes.insert(query_node2);
	left_proprogate.insert(query_node1);
	right_proprogate.insert(query_node2);
	//node1_set.insert(query_node1);
	//node2_set.insert(query_node2);

	//for (auto iter = adjacency_list[query_node1].begin(); iter != adjacency_list[query_node1].end(); iter++)
	//{
	//	left_proprogate.insert(*iter);
	//	reverse_adjacency_in_subgraph_left.insert(std::make_pair(*iter, node1_set));// store the first node
	//}


	//for (auto iter = adjacency_list_reverse[query_node2].begin(); iter != adjacency_list_reverse[query_node2].end(); iter++)
	//{
	//	right_proprogate.insert(*iter);
	//	reverse_adjacency_in_subgraph_right.insert(std::make_pair(*iter, node2_set));// store the first node
	//}

	bool left_skip = false;
	bool right_skip = false;
	int left_max_distance = k;// k - (k / 2);;
	int right_max_distance = k;// (k / 2);
	cur_distance = 0;
	while (true)
	{
		//left_skip = false;
		//right_skip = false;
		cur_distance++;
		if (cur_distance <= left_max_distance && !left_proprogate.empty())
		{
			set<NODE_TYPE> temp_proprogate;
			NODE_TYPE cur_node;
			for (auto iter = left_proprogate.begin(); iter != left_proprogate.end(); iter++)
			{
				cur_node = *iter;

				set<NODE_TYPE> temp_set;
				temp_set.insert(cur_node);
				for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
				{
					if (*iter2 == query_node2)
					{
						result = cur_distance;
						return result;
					}
					//if (right_proprogate.find(*iter2) != right_proprogate.end())// judge whether *iter is a meet node
					//{
					//	meet_nodes.insert(*iter2);
					//	//if (reverse_adjacency_in_subgraph_left.find(*iter2) == reverse_adjacency_in_subgraph_left.end())// not in index subgraph
					//	//{
					//	//	reverse_adjacency_in_subgraph_left.insert(std::make_pair(*iter2, temp_set));
					//	//}
					//	//else {
					//	//	reverse_adjacency_in_subgraph_left.find(*iter2)->second.insert(cur_node);
					//	//}
					//	//	
					//	//continue;
					//}
					set<NODE_TYPE> temp_set2;
					temp_set2.insert(*iter2);
					//if (dag_min_induced_subgraph_reverse.find(*iter2) == dag_min_induced_subgraph_reverse.end())// not in index subgraph
					//{
					//	dag_min_induced_subgraph_reverse.insert(std::make_pair(cur_node, temp_set2));
					//}
					//else {
					//	dag_min_induced_subgraph_reverse.find(cur_node)->second.insert(*iter2);
					//}
					if (left_visited_nodes.find(*iter2) != left_visited_nodes.end())//already visited
					{
						//reverse_adjacency_in_subgraph_left.find(*iter2)->second.insert(cur_node);
						continue;
					}

					left_visited_nodes.insert(*iter2);
					temp_proprogate.insert(*iter2);
				}
			}
			left_proprogate = temp_proprogate;
		}
		else {
			break;//left_skip = true;
		}
	}
	return result;
	//after left bfs, then we do a right bfs to get dag minimum induced subgraph
	//unordered_map<NODE_TYPE, set<NODE_TYPE>> dag_min_induced_subgraph_temp;
	//unordered_map<NODE_TYPE, set<NODE_TYPE>> dag_min_induced_subgraph_reverse_temp;
}


map<NODE_TYPE,DISTANCE_TYPE> bfs_return_spd_pairs(unordered_map<NODE_TYPE, set<NODE_TYPE>> &induced_subgraph, NODE_TYPE  node_num, NODE_TYPE  k, NODE_TYPE query_node1, set<NODE_TYPE>& cut_nodes)
{//store the final induced subgraph into reverse_adjacency_in_subgraph_left
	map<NODE_TYPE, DISTANCE_TYPE> result_pairs;
	int result = -1;
	DISTANCE_TYPE cur_distance = 0;
	set<NODE_TYPE> left_proprogate;
	set<NODE_TYPE> right_proprogate;
	//cur_proprogate.insert(query_node1);
	set<NODE_TYPE> node1_set;
	set<NODE_TYPE> node2_set;
	//set<NODE_TYPE> visited_nodes;

	set<NODE_TYPE> left_visited_nodes;

	//dag_min_induced_subgraph_reverse.insert(std::make_pair(query_node2, node2_set));// query node1's parent is empty
	left_visited_nodes.insert(query_node1);
	left_proprogate.insert(query_node1);
	//node1_set.insert(query_node1);
	//node2_set.insert(query_node2);

	//for (auto iter = adjacency_list[query_node1].begin(); iter != adjacency_list[query_node1].end(); iter++)
	//{
	//	left_proprogate.insert(*iter);
	//	reverse_adjacency_in_subgraph_left.insert(std::make_pair(*iter, node1_set));// store the first node
	//}


	//for (auto iter = adjacency_list_reverse[query_node2].begin(); iter != adjacency_list_reverse[query_node2].end(); iter++)
	//{
	//	right_proprogate.insert(*iter);
	//	reverse_adjacency_in_subgraph_right.insert(std::make_pair(*iter, node2_set));// store the first node
	//}

	bool left_skip = false;
	bool right_skip = false;
	int left_max_distance = k;// k - (k / 2);;
	int right_max_distance = k;// (k / 2);
	cur_distance = 0;
	while (true)
	{
		//left_skip = false;
		//right_skip = false;
		cur_distance++;
		if (cur_distance <= left_max_distance && !left_proprogate.empty())
		{
			set<NODE_TYPE> temp_proprogate;
			NODE_TYPE cur_node;
			for (auto iter = left_proprogate.begin(); iter != left_proprogate.end(); iter++)
			{
				cur_node = *iter;

				set<NODE_TYPE> temp_set;
				temp_set.insert(cur_node);
				auto iter_next = induced_subgraph.find(cur_node);
				if (iter_next == induced_subgraph.end())
				{
					continue;
				}
				for (auto iter2 = iter_next->second.begin(); iter2 != iter_next->second.end(); iter2++)
				{
					if (cut_nodes.find(*iter2) != cut_nodes.end())
					{//insert the node distance pair
						if (result_pairs.find(*iter2) == result_pairs.end())//first time to insert
						{
							result_pairs.insert(std::make_pair(*iter2, cur_distance));
						}
					}
					//if (right_proprogate.find(*iter2) != right_proprogate.end())// judge whether *iter is a meet node
					//{
					//	meet_nodes.insert(*iter2);
					//	//if (reverse_adjacency_in_subgraph_left.find(*iter2) == reverse_adjacency_in_subgraph_left.end())// not in index subgraph
					//	//{
					//	//	reverse_adjacency_in_subgraph_left.insert(std::make_pair(*iter2, temp_set));
					//	//}
					//	//else {
					//	//	reverse_adjacency_in_subgraph_left.find(*iter2)->second.insert(cur_node);
					//	//}
					//	//	
					//	//continue;
					//}
					set<NODE_TYPE> temp_set2;
					temp_set2.insert(*iter2);
					//if (dag_min_induced_subgraph_reverse.find(*iter2) == dag_min_induced_subgraph_reverse.end())// not in index subgraph
					//{
					//	dag_min_induced_subgraph_reverse.insert(std::make_pair(cur_node, temp_set2));
					//}
					//else {
					//	dag_min_induced_subgraph_reverse.find(cur_node)->second.insert(*iter2);
					//}
					if (left_visited_nodes.find(*iter2) != left_visited_nodes.end())//already visited
					{
						//reverse_adjacency_in_subgraph_left.find(*iter2)->second.insert(cur_node);
						continue;
					}

					left_visited_nodes.insert(*iter2);
					temp_proprogate.insert(*iter2);
				}
			}
			left_proprogate = temp_proprogate;
		}
		else {
			break;//left_skip = true;
		}
	}
	return result_pairs;
	//after left bfs, then we do a right bfs to get dag minimum induced subgraph
	//unordered_map<NODE_TYPE, set<NODE_TYPE>> dag_min_induced_subgraph_temp;
	//unordered_map<NODE_TYPE, set<NODE_TYPE>> dag_min_induced_subgraph_reverse_temp;
}



void write_random_query_edges_withspd(string dataset, int k, NODE_TYPE query_node1, NODE_TYPE query_node2, NODE_TYPE edgeIndex, DISTANCE_TYPE spd)
{
	ofstream temp;
	if (edgeIndex == 1)
	{
		temp.open(dataset + "_" + to_string(k)+"_randomquery");
	}
	else
	{
		temp.open(dataset + "_" + to_string(k)+"_randomquery", ios::app);
	}
	if (!temp)
	{
		// cout << "random querys can't open" << endl;
		abort();
	}
	temp << query_node1 << "\t" << query_node2 << "\t" << spd << endl;
	temp.close();
}




void write_random_query_edges(string dataset, int k, NODE_TYPE query_node1, NODE_TYPE query_node2, NODE_TYPE edgeIndex)
{
	ofstream temp;
	if (edgeIndex == 1)
	{
		temp.open(dataset + "_" + to_string(k) + "_randomquery");
	}
	else
	{
		temp.open(dataset + "_" + to_string(k) + "_randomquery", ios::app);
	}
	if (!temp)
	{
		// cout << "random querys can't open" << endl;
		abort();
	}
	temp << query_node1 << "\t" << query_node2 << endl;
	temp.close();
}


void random_query_pair(NODE_TYPE & query_node1, NODE_TYPE &query_node2, NODE_TYPE nodenum)
{
	std::random_device rd;
	std::uniform_int_distribution<int> dist(1, nodenum);
	query_node1 = dist(rd);
	query_node2 = dist(rd);
	while (query_node2 == query_node1)
	{
		query_node2 = dist(rd);
	}
}

bool spd_index::judge_reachability(NODE_TYPE node1, NODE_TYPE node2, DISTANCE_TYPE distance)//node1 must be in hot points
{
	bool result = false;
	auto iter1 = h_nh_index.find(node1);
	if (iter1 != h_nh_index.end())//  non-hot index for node1 and can reach node2 in distance contraint
	{
		auto iter2 = iter1->second.distance_map.find(node2);
		if (iter2 != iter1->second.distance_map.end() && iter2->second <= distance)// node2 is reachability from node1 under distance contraint
		{
			return true;
		}
	}
	auto iter3 = h_h_index.find(node1);
	if (iter3 == h_h_index.end())// no hot-hot index for node1
	{
		return false;
	}
	else {
		auto iter2 = iter3->second.distance_map.find(node2);
		if (iter2 != iter3->second.distance_map.end() && iter2->second <= distance)// node2 is reachability from node1 under distance contraint
		{
			return true;
		}

		for (auto iter4 = iter3->second.distance_map.begin(); iter4 != iter3->second.distance_map.end(); iter4++)
		{
			if (iter4->second <= distance) //distance <= remaining distance
			{
				result = result || judge_reachability(iter4->first, node2, distance - iter4->second);
			}
		}
	}
	return result;
}

void spd_index::construct_spd_index(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  k)
{
	if (h_h_index.size() != 0 || h_nh_index.size() != 0)
	{
		// cout << "already has an index, clear it first" << endl;
		return;
	}
	// cout << " k is " << k << " in contruct index " << endl;
	// use bfs to find some paths between node A and B
	for (auto start = hot_points.begin(); start != hot_points.end(); start++)
	{
		//spd_distance_map temp_set;
		//spd_distance_map temp;
		//temp = single_direction_baseline(adjacency_list,node_num,*start,*end,k);
		single_direction_baseline_stop_at_hotpoints_bfs_spd(adjacency_list, node_num, *start, k, hot_points, *this);
		//temp_set.add_map(temp);

		//h_h_index.insert(std::make_pair(*start, temp_set));
	}
	return;
	// finish index construct

}

paths find_k_paths_between_two_nodes_spd(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE k, NODE_TYPE query_node1, NODE_TYPE query_node2, spd_index& index)
{
	paths result;
	int min_stop_distance = int(k);
	int right_stop_distance = int(k);

	// some structures to store unhot paths
	set<NODE_TYPE> stop_points(index.hot_points);
	stop_points.insert(query_node2);
	map<NODE_TYPE, paths> left_dfs;
	map<NODE_TYPE, paths> right_dfs;
	vector<NODE_TYPE> c_path;
	set<NODE_TYPE> c_path_set;
	dfs(adjacency_list_double, adjacency_list, left_dfs, query_node1, c_path, stop_points, k, 0, min_stop_distance, c_path_set);

	paths paths_without_hot;
	auto iter_left = left_dfs.find(query_node2);
	if (iter_left != left_dfs.end())
	{
		paths_without_hot.add_paths(iter_left->second);
		//paths_without_hot.output();
	}
	result.add_paths(paths_without_hot);
	for (auto iter1 = left_dfs.begin(); iter1 != left_dfs.end(); iter1++)
	{
		NODE_TYPE cur_node = iter1->first;
		if (cur_node == query_node2)
		{
			continue;
		}
		if (index.judge_reachability(cur_node, query_node2, k - iter1->second.get_min_distance()))//can reach
		{
			paths temp, temp_result;
			set<NODE_TYPE> cur_set;
			vector<NODE_TYPE> cur_path;
			dfs_find_k_paths(adjacency_list, cur_node, cur_node, query_node2, k - iter1->second.get_min_distance(), temp, cur_set, 0, cur_path);
			temp_result = iter1->second.join(temp);
			result.add_paths(temp_result);
		}
		else {
			continue;
		}
	}
	result.drop_path_length_more_than_k(k);
	result.drop_path_with_repeat_node();
	result.sort_by_string_order();
	return result;
}


paths start_paths::get_paths()
{
	paths result;
	for(auto iter = t_path.begin(); iter != t_path.end(); iter++)
	{
		paths temp;
		temp = iter->second.get_paths();
		result.add_paths(temp);
	}
	return result;
}



void extractEdges(char str[], NODE_TYPE  &x, NODE_TYPE  &y)
{
	//char* end;
	char temp[256];
	NODE_TYPE  tempIndex = 0;
	NODE_TYPE  i = 0;
	for (; '0' <= str[i] && str[i] <= '9'; i++)
	{
		temp[tempIndex++] = str[i];
	}
	temp[tempIndex] = '\0';
	x = atoi(temp);
	tempIndex = 0;
	while (str[i + 1] != '\0' && str[i + 1] != ' ' && str[i + 1] != '\t')
	{
		temp[tempIndex++] = str[++i];
	}
	temp[tempIndex] = '\0';
	y = atoi(temp);
	return;
}

void extractEdgesSkipNodes(char str[], NODE_TYPE  &x, NODE_TYPE  &y,int skip_node)
{
	//char* end;
	char temp[256];
	NODE_TYPE  tempIndex = 0;
	NODE_TYPE  i = 0;
	for (int j = 0; j < skip_node; j++)
	{
		for (;str[i] != ' ' && str[i] != '\t' && str[i] != '\0';i++)
		{
			//i++;
		}
		i++;
	}
	for (; '0' <= str[i] && str[i] <= '9'; i++)
	{
		temp[tempIndex++] = str[i];
	}
	temp[tempIndex] = '\0';
	x = atoi(temp);
	tempIndex = 0;
	while (str[i+1] != '\0' && str[i+1] != ' ' && str[i+1] != '\t')
	{
		temp[tempIndex++] = str[++i];
	}
	temp[tempIndex] = '\0';
	y = atoi(temp);
	return;
}

void extractAliEdges(char str[], NODE_TYPE &x, NODE_TYPE &y)
{
	char temp[256];
	NODE_TYPE  tempIndex = 0;
	NODE_TYPE  i = 0;
	while (str[i++] != ',') {
		;
	}
	for (; '0' <= str[i] && str[i] <= '9'; i++)
	{
		temp[tempIndex++] = str[i];
	}
	temp[tempIndex] = '\0';
	istringstream ss(temp);
	ss >> x;
	while (str[i++] != ',') {
		;
	}
	while (str[i++] != ',') {
		;
	}
	tempIndex = 0;
	while (str[++i] != ',')
	{
		temp[tempIndex++] = str[i];
	}
	temp[tempIndex] = '\0';
	istringstream ss1(temp);
	ss1 >> y;
	return;
}


void load_ali_undirected_data(vector<NODE_TYPE>* adjacency_list_double, vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse)// if the graph has an undirected edge (a,b), there should be (a,b) and (b,a) in the record file
{

	int table_range[6] = { 0,7,11,11,11,7 };//table from 1 to 5 and range from 0 to table range[i]
	int table = 1;
	NODE_TYPE x, y;
	string table_name;
	ifstream inEdges;
	char buffer[BUFFER_LENTH];
	//string x,y;
	for (table = 0; table <= 19; table++)
	{
		table_name = "./data/dwd_rsm_dmt_graph_data20171228_01_dd_idx/dwd_rsm_dmt_graph_data20171228_01_dd_idx_" + to_string(table);
		// cout << table_name << " is our talbe name" << endl;
		//for each table, input all the key-value pairs. Store and output all ordered key-value pairs.
		inEdges.open(table_name);
		if (!inEdges.is_open())
		{
			// cout << "Error opening file : " << table_name << endl;
			continue;
		}
		while (!inEdges.eof())
		{
			inEdges.getline(buffer, BUFFER_LENTH);
			extractEdges(buffer, x, y);
			adjacency_list_double[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list
			//adjacency_list_reverse[y].push_back(x);
		}
		inEdges.close();
	}

	for (table = 0; table <= 19; table++)
	{
		table_name = "./data/dwd_rsm_dmt_graph_data20171228_04_two_days_dd_idx/dwd_rsm_dmt_graph_data20171228_04_two_days_dd_idx_" + to_string(table);
		// cout << table_name << " is our talbe name" << endl;
		//for each table, input all the key-value pairs. Store and output all ordered key-value pairs.
		inEdges.open(table_name);
		if (!inEdges.is_open())
		{
			// cout << "Error opening file : " << table_name << endl;
			continue;
		}
		while (!inEdges.eof())
		{
			inEdges.getline(buffer, BUFFER_LENTH);
			extractEdges(buffer, x, y);
			adjacency_list[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list
			adjacency_list_reverse[y].push_back(x);					  //adjacency_list_reverse[y].push_back(x);
		}
		inEdges.close();
	}


	return;
}

void load_ali_data(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse)
{

	int table_range[6] = { 0,7,11,11,11,7 };//table from 1 to 5 and range from 0 to table range[i]
	int table = 1;
	NODE_TYPE x, y;
	string table_name;
	ifstream inEdges;
	char buffer[BUFFER_LENTH];
	//string x,y;
	for (table = 1; table <= 5; table++)
	{
		for (int cur_index = 0; cur_index <= table_range[table]; cur_index++)
		{

			table_name = "table" + to_string(table) + "/table" + to_string(table) + "_" + to_string(cur_index) + "out.txt";
			// cout << table_name << " is our talbe name" << endl;
			//for each table, input all the key-value pairs. Store and output all ordered key-value pairs.
			inEdges.open(table_name);
			if (!inEdges.is_open())
			{
				// cout << "Error opening file : " << table_name << endl;
				continue;
			}
			while (!inEdges.eof())
			{
				inEdges.getline(buffer, BUFFER_LENTH);
				extractEdges(buffer, x, y);
				adjacency_list[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list
				adjacency_list_reverse[y].push_back(x);
			}
			inEdges.close();
		}
	}
	return ;
}

void load_graph_data(string filename, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse)
{
	char buffer[BUFFER_LENTH];
	ifstream inEdges(filename);
	if (!inEdges.is_open())
	{
		// cout << "Error opening file" << endl;
		exit(1);
	}
	NODE_TYPE  x, y;
	int i = 0;
	bool roadUSA = false;
	bool pound = false;
	while (!inEdges.eof())
	{
		i++;
		inEdges.getline(buffer, BUFFER_LENTH);
		if (i == 1 && buffer[0] == 'c')
		{
			roadUSA = true;
		}
		else if (i == 1 && (buffer[0] == '#' || (buffer[0] >= '0' && buffer[0] <= '9')))
		{
			pound = true;
		}
		if (roadUSA)
		{
			if (buffer[0] == 'c' || buffer[0] == 'p')
			{
				continue;
			}

			extractEdgesSkipNodes(buffer, x, y, 1);
			if (x == 0 || y == 0) continue;
			if (x == y)
			{
				continue;//remove self-loop
			}
			adjacency_list[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list

			adjacency_list_reverse[y].push_back(x);
		}
		else if(pound)
		{
			if (buffer[0] == '#')
			{
				continue;
			}

			extractEdges(buffer, x, y);
			if (x == 0 || y == 0) continue;
			if (x == y)
			{
				continue;//remove self-loop
			}
			adjacency_list[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list

			adjacency_list_reverse[y].push_back(x);
			//adjacency_list[y].push_back(x);
			//adjacency_list_reverse[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list
		}
		
		//adjacency_list[y].push_back(x);
		//adjacency_list_reverse[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list


	}
	return;
}




void load_graph_data(char* filename, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse)
{
	char buffer[BUFFER_LENTH];
	ifstream inEdges(filename);
	if (!inEdges.is_open())
	{
		// cout << "Error opening file" << endl;
		exit(1);
	}
	NODE_TYPE  x, y;
	int i = 0;
	while (!inEdges.eof())
	{
		i++;
		inEdges.getline(buffer, BUFFER_LENTH);
		if (i <= 10)
		{
			continue;
		}
		extractEdges(buffer, x, y);
		if (x == y)
		{
			continue;
		}
		adjacency_list[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list
		adjacency_list_reverse[y].push_back(x);
		//adjacency_list[y].push_back(x);
	}
	return;
}


void load_undirected_graph_data(char* filename, vector<NODE_TYPE >* adjacency_list_double)
{
	char buffer[BUFFER_LENTH];
	ifstream inEdges(filename);
	if (!inEdges.is_open())
	{
		// cout << "Error opening file" << endl;
		exit(1);
	}
	NODE_TYPE  x, y;
	int i = 0;
	while (!inEdges.eof())
	{
		i++;
		inEdges.getline(buffer, BUFFER_LENTH);
		if (i <= 10)
		{
			continue;
		}
		extractEdges(buffer, x, y);
		adjacency_list_double[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list
		//adjacency_list_reverse[y].push_back(x);
		//adjacency_list[y].push_back(x);
	}
	return;
}

void output_adjacency_list(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num)
{
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		if (adjacency_list[i].size() >= 1)
		{
			// cout << i << " ";
			for (auto iter = adjacency_list[i].begin(); iter != adjacency_list[i].end(); iter++)
			{
				// cout << (*iter) << " ";
			}
			// cout << endl;
		}
	}
}

void remove_and_output_meetpath(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, map<NODE_TYPE, vector<NODE_TYPE > >* partents,
	map<NODE_TYPE, vector<NODE_TYPE > >* parents_meetup)
{

}


void double_direction_unordered(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k)
{// use unordered_map to finish double_direction algorithm
	unordered_map<NODE_TYPE, vector<string> > first_path_map;
	unordered_map<NODE_TYPE, vector<string> > second_path_map;
	ofstream out_path;
	out_path.open("double_out_unordered.txt");
	//set<pair<NODE_TYPE ,NODE_TYPE > >* parents = new set<pair<NODE_TYPE ,NODE_TYPE > >[node_num+1];
	NODE_TYPE path_count = 0;
	bool* is_board_nodes = new bool[node_num + 1];
	for (int i = 0; i <= node_num; i++)
	{
		is_board_nodes[i] = false;
	}
	unordered_map<NODE_TYPE, NODE_TYPE> path_count_map;
	unordered_map<NODE_TYPE, NODE_TYPE> path_count_map_second;
	unordered_map<NODE_TYPE, vector<NODE_TYPE > >* parents = new unordered_map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];//record distance and last nodes
																									   //map<NODE_TYPE, vector<NODE_TYPE > >* parents_meetup = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];//record distance and last nodes
	set<NODE_TYPE > cur_propagate;
	NODE_TYPE  cur_distance = 1;
	cur_propagate.insert(query_node1);
	while (cur_distance <= int(k / 2))
	{
		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			//// cout << "cur_node is " << cur_node << endl;
			for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes
			 //// cout << "iter2 " << *iter2 << endl;

				if (parents[*iter2].find(cur_distance) == parents[*iter2].end()) {
					//if there is no vector for this distance, pair is distance, node 
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
					//parents_meetup[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
													  //iter3 = parents_meetup[*iter2].find(cur_distance);
													  //iter3->second.push_back(cur_node);//insert node distance value
				}
				//		auto iter4 = parents[*iter2].find(cur_distance);
				if (cur_distance == int(k / 2))//record all the board nodes
				{
					is_board_nodes[*iter2] = true;
				}
				temp_propagate.insert(*iter2);//insert all the new propagate nodes
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}

	cur_propagate.clear();
	cur_propagate.insert(query_node2);//search reversely
	cur_distance = 1;
	unordered_map<NODE_TYPE, vector<NODE_TYPE > >* parents_query2 = new unordered_map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];//record distance and last nodes
	stack<pair<NODE_TYPE, NODE_TYPE > > short_path;
	if (!parents[query_node2].empty())
	{
		//output all the short path and remove all the nodes on these path
		for (auto iter = parents[query_node2].begin(); iter != parents[query_node2].end(); iter++)
		{
			for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
			{
				short_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter->first));
				//parents_meetup[query_node2].erase(iter->first);

			}
		}
	}
	vector<NODE_TYPE> output_path;
	output_path.push_back(query_node2);
	if (first_path_map.find(query_node2) == first_path_map.end())
	{
		vector<string> a;
		first_path_map.insert(pair<NODE_TYPE, vector<string> >(query_node2, a));
		a.push_back("");
		second_path_map.insert(pair<NODE_TYPE, vector<string> >(query_node2, a));
	}
	string temp_path_map("");

	// cout << "short path is " << endl;
	while (!short_path.empty())
	{
		NODE_TYPE  temp_node = short_path.top().first;
		NODE_TYPE  distance = short_path.top().second;
		if (temp_node == query_node1)
		{
			//output_path
			temp_path_map.clear();
			for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
			{
				// cout << *iter4 << " ";
				//out_path << *iter4 << " ";
				temp_path_map += to_string(*iter4) + " ";
			}
			// cout << query_node1 << endl;
			//out_path << query_node1 << endl;
			temp_path_map += to_string(query_node1);
			first_path_map.find(query_node2)->second.push_back(temp_path_map);
			path_count++;
			//output_path.pop_back();
			short_path.pop();
			continue;
		}
		else if (temp_node == output_path.back())
		{
			output_path.pop_back();
			short_path.pop();
			continue;
		}
		else {
			output_path.push_back(temp_node);
		}
		//path.pop();
		/*if (distance == k)
		{
		// cout << endl;
		}*/
		//// cout << temp_node<<":" << distance << " ";
		distance--;
		auto iter = parents[temp_node].find(distance);
		/*if (distance == 0 || iter == parents[temp_node].end())
		{
		// cout << " end of path" << endl;
		continue;
		}*/
		auto iter2 = iter->second;
		for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
		{
			short_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
		}
	}
	set<NODE_TYPE > another_path;
	is_board_nodes[query_node2] = false;// because query_node2 is outputed in the short path section
	is_board_nodes[query_node1] = false;
	while (cur_distance <= k - int(k / 2))// meet up reversely 
	{
		//start to output all the meetup paths

		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			//// cout << "cur_node is " << cur_node << endl;
			for (auto iter2 = adjacency_list_reverse[cur_node].begin(); iter2 != adjacency_list_reverse[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes

				if (parents_query2[*iter2].find(cur_distance) == parents_query2[*iter2].end()) {
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents_query2[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents_query2[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
				}
				//		auto iter4 = parents[*iter2].find(cur_distance);

				temp_propagate.insert(*iter2);//insert all the new propagate nodes
											  //judge whether there are meet-up paths at cur_node
				if (is_board_nodes[*iter2] == true)
				{
					if (first_path_map.find(*iter2) == first_path_map.end())
					{
						vector<string> a;
						first_path_map.insert(pair<NODE_TYPE, vector<string> >(*iter2, a));
						second_path_map.insert(pair<NODE_TYPE, vector<string> >(*iter2, a));
					}
					string temp_path_map("");
					if (*iter2 == 12647)
					{
						// cout << 12647 << endl;
					}
					//output the first half path
					is_board_nodes[*iter2] = false;
					// cout << "meetup node is " << *iter2 << endl;
					//out_path << "meetup node is " << *iter2 << endl;
					path_count_map.insert(pair<NODE_TYPE, NODE_TYPE>(*iter2, 0));
					path_count_map_second.insert(pair<NODE_TYPE, NODE_TYPE>(*iter2, 0));
					// cout << endl << "first half path: " << endl;
					//out_path << endl << "first half path: " << endl;
					stack<pair<NODE_TYPE, NODE_TYPE > > temp_path;
					while (!temp_path.empty())
					{
						temp_path.pop();
					}
					//output all the short path and remove all the nodes on these path
					auto iter_temp_board = parents[*iter2].find(int(k / 2));
					if (iter_temp_board != parents[*iter2].end())
					{
						for (auto iter_board = iter_temp_board->second.begin(); iter_board != iter_temp_board->second.end(); iter_board++)
						{
							temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter_board, int(k / 2)));
						}
					}

					//temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2,int(k/2)+1 ) ); // distacne node pair
					/*for (auto iter_meet = parents_meetup[cur_node].begin(); iter_meet != parents_meetup[cur_node].end(); iter++)
					{
					for (auto iter2 = iter_meet->second.begin(); iter2 != iter_meet->second.end(); iter2++)
					{
					temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter_meet->first));
					}
					}*/
					//parents_meetup[*iter2].clear();
					auto temp_iter = path_count_map.find(*iter2);
					vector<NODE_TYPE> output_path;
					//output_path.push_back(*iter2);// to be confirm
					while (!temp_path.empty())
					{
						NODE_TYPE  temp_node = temp_path.top().first;
						NODE_TYPE  distance = temp_path.top().second;
						if (temp_node == query_node1)
						{
							//output_path
							temp_path_map.clear();
							//if (!circle_detection(output_path))
							{
								for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
								{
									// cout << *iter4 << " ";
									//out_path << *iter4 << " ";
									temp_path_map += to_string(*iter4) + " ";
								}
								// cout << query_node1 << endl;
								//out_path << query_node1 << endl;
								temp_path_map += to_string(query_node1);
								first_path_map.find(*iter2)->second.push_back(temp_path_map);
								temp_iter->second += 1;
								//output_path.pop_back();
							}
							temp_path.pop();
							continue;
						}
						else if (output_path.size() == 0)
						{
							output_path.push_back(temp_node);
						}
						else if (temp_node == output_path.back())
						{
							output_path.pop_back();
							temp_path.pop();
							continue;
						}
						else {
							output_path.push_back(temp_node);
						}
						//temp_path.pop();
						//// cout << temp_node << " ";
						distance--;
						if (parents[temp_node].find(distance) == parents[temp_node].end())
						{
							while (!temp_path.empty())
							{
								temp_path.pop();
							}
							break;
						}
						auto iter_parent = parents[temp_node].find(distance);
						//auto iter3 = parents_meetup[temp_node].find(distance);
						/*if (distance == 0 || iter == parents[temp_node].end())
						{
						// cout << " end of path" << endl;
						continue;
						}*/
						auto iter2 = iter_parent->second;
						//iter3->second.clear();
						//parents_meetup[temp_node].erase(distance);
						for (auto iter2 = iter_parent->second.begin(); iter2 != iter_parent->second.end(); iter2++)
						{
							temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
						}
					}

					// output another half path
					another_path.insert(*iter2);

				}
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}

	// cout << "another half path" << endl;
	//out_path << "another half path" << endl;
	stack<pair<NODE_TYPE, NODE_TYPE > > temp_path;
	for (auto iter_set = another_path.begin(); iter_set != another_path.end(); iter_set++)
	{
		NODE_TYPE  cur_node = *iter_set;
		// cout << "meet node is " << cur_node << endl;
		//out_path << "meet node is " << cur_node << endl;
		if (!parents_query2[cur_node].empty())
		{
			//output all the short path and remove all the nodes on these path
			for (auto iter = parents_query2[cur_node].begin(); iter != parents_query2[cur_node].end(); iter++)
			{
				for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
				{
					temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter->first));
				}
			}
		}
		auto temp_iter2 = path_count_map_second.find(cur_node);
		vector<NODE_TYPE> output_path;
		output_path.push_back(cur_node);
		string temp_path_map("");
		while (!temp_path.empty())
		{
			NODE_TYPE  temp_node = temp_path.top().first;
			NODE_TYPE  distance = temp_path.top().second;
			//temp_path.pop();
			//// cout << temp_node << " ";
			if (temp_node == query_node2)
			{
				temp_path_map.clear();
				//output_path
				for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
				{
					// cout << *iter4 << " ";
					//out_path << *iter4 << " ";
					temp_path_map += to_string(*iter4) + " ";
				}
				// cout << query_node2 << endl;
				//out_path << query_node2 << endl;
				temp_path_map += to_string(query_node2);
				second_path_map.find(cur_node)->second.push_back(temp_path_map);
				temp_iter2->second += 1;
				//output_path.pop_back();
				temp_path.pop();
				continue;
			}
			else if (temp_node == output_path.back())
			{
				output_path.pop_back();
				temp_path.pop();
				continue;
			}
			else {
				output_path.push_back(temp_node);
			}
			distance--;
			auto iter = parents_query2[temp_node].find(distance);
			//auto iter3 = parents_meetup[temp_node].find(distance);
			/*if (distance == 0 || iter == parents_query2[temp_node].end())
			{
			// cout << " end of path" << endl;
			continue;
			}*/
			auto iter2 = iter->second;
			//iter3->second.clear();
			for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
			{
				temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
			}
		}
	}
	for (auto iter_count = path_count_map.begin(); iter_count != path_count_map.end(); iter_count++)
	{
		auto iter_count2 = path_count_map_second.find(iter_count->first);
		//// cout << iter_count->first << " is meet point" << iter_count2->second << ":" << iter_count->second << " in total" << endl;
		path_count += iter_count2->second * iter_count->second;
	}

	// cout << "all the path" << endl;
	//out_path << "all the path" << endl;
	for (auto iter = first_path_map.begin(); iter != first_path_map.end(); iter++)
	{
		auto iter2 = second_path_map.find(iter->first);
		for (auto iter3 = iter->second.begin(); iter3 != iter->second.end(); iter3++)
		{
			for (auto iter4 = iter2->second.begin(); iter4 != iter2->second.end(); iter4++)
			{
				if (!circle_detection_string(reverse_path(*iter3) + *iter4))
				{
					// cout << reverse_path(*iter3) << *iter4 << endl;
					out_path << reverse_path(*iter3) << *iter4 << " " << endl;
				}
			}
		}
	}


	out_path.close();
	// cout << path_count << " path in double direction method" << endl;
	return;

}




void double_direction_baseline(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k)
{  
	map<NODE_TYPE, vector<string> > first_path_map;
	map<NODE_TYPE, vector<string> > second_path_map;
	ofstream out_path;
	out_path.open("double_bfs_result",ios::app);
	NODE_TYPE path_count = 0;
	bool* is_board_nodes = new bool[node_num+1];
	for (int i = 0; i <= node_num; i++)
	{
		is_board_nodes[i] = false;
	}
	map<NODE_TYPE, NODE_TYPE> path_count_map;
	map<NODE_TYPE, NODE_TYPE> path_count_map_second;
	map<NODE_TYPE, vector<NODE_TYPE > >* parents = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];//record distance and last nodes
	set<NODE_TYPE > cur_propagate;
	NODE_TYPE  cur_distance = 1;
	cur_propagate.insert(query_node1);
	while (cur_distance <= int(k / 2))
	{
		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node =  *iter;
			//// cout << "cur_node is " << cur_node << endl;
			for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes
			 //// cout << "iter2 " << *iter2 << endl;

				if (parents[*iter2].find(cur_distance) == parents[*iter2].end()) {
					//if there is no vector for this distance, pair is distance, node 
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
					//parents_meetup[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
					//iter3 = parents_meetup[*iter2].find(cur_distance);
					//iter3->second.push_back(cur_node);//insert node distance value
				}
				//		auto iter4 = parents[*iter2].find(cur_distance);
				if (cur_distance == int(k / 2))//record all the board nodes
				{
					is_board_nodes[*iter2] = true;
				}
				temp_propagate.insert(*iter2);//insert all the new propagate nodes
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}

	cur_propagate.clear();
	cur_propagate.insert(query_node2);//search reversely
	cur_distance = 1;
	map<NODE_TYPE, vector<NODE_TYPE > >* parents_query2 = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];//record distance and last nodes
	stack<pair<NODE_TYPE, NODE_TYPE > > short_path;
	if (!parents[query_node2].empty())
	{
		//output all the short path and remove all the nodes on these path
		for (auto iter = parents[query_node2].begin(); iter != parents[query_node2].end(); iter++)
		{
			for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
			{
				short_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter->first));
				//parents_meetup[query_node2].erase(iter->first);

			}
		}
	}
	vector<NODE_TYPE> output_path;
	output_path.push_back(query_node2);
	if (first_path_map.find(query_node2)  == first_path_map.end())
	{
		vector<string> a;
		first_path_map.insert(pair<NODE_TYPE, vector<string> >(query_node2,a));
		a.push_back("");
		second_path_map.insert(pair<NODE_TYPE, vector<string> >(query_node2, a));
	}
	string temp_path_map("");

	//// cout << "short path is " << endl;
	while (!short_path.empty())
	{
		NODE_TYPE  temp_node = short_path.top().first;
		NODE_TYPE  distance = short_path.top().second;
		if (temp_node == query_node1)
		{
			//output_path
			temp_path_map.clear();
			for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
			{
				//// cout << *iter4 << " ";
				//out_path << *iter4 << " ";
				temp_path_map += to_string(*iter4) + " ";
			}
			//// cout << query_node1 << endl;
			//out_path << query_node1 << endl;
			temp_path_map += to_string(query_node1);
			first_path_map.find(query_node2)->second.push_back(temp_path_map);
			path_count++;
			//output_path.pop_back();
			short_path.pop();
			continue;
		}
		else if (temp_node == output_path.back())
		{
			output_path.pop_back();
			short_path.pop();
			continue;
		}
		else {
			output_path.push_back(temp_node);
		}
		//path.pop();
		/*if (distance == k)
		{
		// cout << endl;
		}*/
		//// cout << temp_node<<":" << distance << " ";
		distance--;
		auto iter = parents[temp_node].find(distance);
		/*if (distance == 0 || iter == parents[temp_node].end())
		{
		// cout << " end of path" << endl;
		continue;
		}*/
		auto iter2 = iter->second;
		for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
		{
			short_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
		}
	}
	set<NODE_TYPE > another_path;
	is_board_nodes[query_node2] = false;// because query_node2 is outputed in the short path section
	is_board_nodes[query_node1] = false;
	while (cur_distance <= k - int(k / 2))// meet up reversely 
	{
		//start to output all the meetup paths

		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
 			NODE_TYPE  cur_node = *iter;
			//// cout << "cur_node is " << cur_node << endl;
			for (auto iter2 = adjacency_list_reverse[cur_node].begin(); iter2 != adjacency_list_reverse[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes

				if (parents_query2[*iter2].find(cur_distance) == parents_query2[*iter2].end()) {
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents_query2[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents_query2[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
				}
				//		auto iter4 = parents[*iter2].find(cur_distance);

				temp_propagate.insert(*iter2);//insert all the new propagate nodes
				//judge whether there are meet-up paths at cur_node
				if(is_board_nodes[*iter2] == true)
				{
					if (first_path_map.find(*iter2) == first_path_map.end())
					{
						vector<string> a;
						first_path_map.insert(pair<NODE_TYPE, vector<string> >(*iter2, a));
						second_path_map.insert(pair<NODE_TYPE, vector<string> >(*iter2, a));
					}
					string temp_path_map("");
					//output the first half path
					is_board_nodes[*iter2] = false;
					//// cout << "meetup node is " << *iter2 << endl;
					//out_path << "meetup node is " << *iter2 << endl;
					path_count_map.insert(pair<NODE_TYPE, NODE_TYPE>(*iter2,0));
					path_count_map_second.insert(pair<NODE_TYPE, NODE_TYPE>(*iter2, 0));
					//// cout << endl << "first half path: " << endl;
					//out_path << endl << "first half path: " << endl;
					stack<pair<NODE_TYPE, NODE_TYPE > > temp_path;
					while (!temp_path.empty())
					{
						temp_path.pop();
					}
					//output all the short path and remove all the nodes on these path
					auto iter_temp_board = parents[*iter2].find(int(k/2));
					if (iter_temp_board != parents[*iter2].end())
					{
						for (auto iter_board = iter_temp_board->second.begin(); iter_board != iter_temp_board->second.end(); iter_board++)
						{
							temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter_board, int(k / 2)));
						}
					}

					//temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2,int(k/2)+1 ) ); // distacne node pair
					/*for (auto iter_meet = parents_meetup[cur_node].begin(); iter_meet != parents_meetup[cur_node].end(); iter++)
					{
						for (auto iter2 = iter_meet->second.begin(); iter2 != iter_meet->second.end(); iter2++)
						{
							temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter_meet->first));
						}
					}*/
					//parents_meetup[*iter2].clear();
					auto temp_iter = path_count_map.find(*iter2);
					vector<NODE_TYPE> output_path;
					//output_path.push_back(*iter2);// to be confirm
					while (!temp_path.empty())
					{
						NODE_TYPE  temp_node = temp_path.top().first;
						NODE_TYPE  distance = temp_path.top().second;
						if (temp_node == query_node1 )
						{
							//output_path
							temp_path_map.clear();
							//if (!circle_detection(output_path))
							{
								for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
								{
									//// cout << *iter4 << " ";
									//out_path << *iter4 << " ";
									temp_path_map += to_string(*iter4) + " ";
								}
								//// cout << query_node1 << endl;
								//out_path << query_node1 << endl;
								temp_path_map += to_string(query_node1);
								first_path_map.find(*iter2)->second.push_back(temp_path_map);
								temp_iter->second += 1;
								//output_path.pop_back();
							}
							temp_path.pop();
							continue;
						}
						else if (output_path.size() == 0)
						{
							output_path.push_back(temp_node);
						}
						else if (temp_node == output_path.back())
						{
							output_path.pop_back();
							temp_path.pop();
							continue;
						}
						else {
							output_path.push_back(temp_node);
						}
						//temp_path.pop();
						//// cout << temp_node << " ";
						distance--;
						if (parents[temp_node].find(distance) == parents[temp_node].end())
						{
							while (!temp_path.empty())
							{
								temp_path.pop();
							}
							break;
						}
						auto iter_parent = parents[temp_node].find(distance);
						//auto iter3 = parents_meetup[temp_node].find(distance);
						/*if (distance == 0 || iter == parents[temp_node].end())
						{
							// cout << " end of path" << endl;
							continue;
						}*/
						auto iter2 = iter_parent->second;
						//iter3->second.clear();
						//parents_meetup[temp_node].erase(distance);
						for (auto iter2 = iter_parent->second.begin(); iter2 != iter_parent->second.end(); iter2++)
						{
							temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
						}
					}

					// output another half path
					another_path.insert(*iter2);

				}
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}

	//// cout << "another half path" << endl;
	//out_path << "another half path" << endl;
	stack<pair<NODE_TYPE, NODE_TYPE > > temp_path;
	for (auto iter_set = another_path.begin(); iter_set != another_path.end(); iter_set++)
	{
		NODE_TYPE  cur_node = *iter_set;
		//// cout << "meet node is " << cur_node << endl;
		//out_path << "meet node is " << cur_node << endl;
		if (!parents_query2[cur_node].empty())
		{
			//output all the short path and remove all the nodes on these path
			for (auto iter = parents_query2[cur_node].begin(); iter != parents_query2[cur_node].end(); iter++)
			{
				for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
				{
					temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter->first));
				}
			}
		}
		auto temp_iter2 = path_count_map_second.find(cur_node);
		vector<NODE_TYPE> output_path;
		output_path.push_back(cur_node);
		string temp_path_map("");
		while (!temp_path.empty())
		{
			NODE_TYPE  temp_node = temp_path.top().first;
			NODE_TYPE  distance = temp_path.top().second;
			//temp_path.pop();
			//// cout << temp_node << " ";
			if (temp_node == query_node2)
			{
				temp_path_map.clear();
				//output_path
				for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
				{
					//// cout << *iter4 << " ";
					//out_path << *iter4 << " ";
					temp_path_map += to_string(*iter4) + " ";
				}
				//// cout << query_node2 << endl;
				//out_path << query_node2 << endl;
				temp_path_map += to_string(query_node2);
				second_path_map.find(cur_node)->second.push_back(temp_path_map);
				temp_iter2->second += 1;
				//output_path.pop_back();
				temp_path.pop();
				continue;
			}
			else if (temp_node == output_path.back())
			{
				output_path.pop_back();
				temp_path.pop();
				continue;
			}
			else {
				output_path.push_back(temp_node);
			}
			distance--;
			auto iter = parents_query2[temp_node].find(distance);
			//auto iter3 = parents_meetup[temp_node].find(distance);
			/*if (distance == 0 || iter == parents_query2[temp_node].end())
			{
				// cout << " end of path" << endl;
				continue;
			}*/
			auto iter2 = iter->second;
			//iter3->second.clear();
			for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
			{
				temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
			}
		}
	}
	for (auto iter_count = path_count_map.begin(); iter_count != path_count_map.end(); iter_count++)
	{
		auto iter_count2 = path_count_map_second.find(iter_count->first);
		//// cout << iter_count->first << " is meet point" << iter_count2->second << ":" << iter_count->second << " in total" << endl;
		path_count += iter_count2->second * iter_count->second;
	}

	//// cout << "all the path" << endl;
	//out_path << "all the path" << endl;
	for (auto iter = first_path_map.begin(); iter != first_path_map.end(); iter++)
	{
		auto iter2 = second_path_map.find(iter->first);
		for (auto iter3 = iter->second.begin(); iter3 != iter->second.end(); iter3++)
		{
			for (auto iter4 = iter2->second.begin(); iter4 != iter2->second.end(); iter4++)
			{
				if (!circle_detection_string(reverse_path(*iter3)  + *iter4))
				{
					//// cout << reverse_path(*iter3) << *iter4 << endl;
					//vector<NODE_TYPE> temp_result;
					//vector<NODE_TYPE> temp2 = reverse_path(*iter3);
					//temp_result.insert(temp_result.end()), 
					//result.push_back(reverse_path(*iter3) + (*iter4));
					out_path << reverse_path(*iter3)<< *iter4 << " " << endl;
				}
			}
		}
	}

	
	out_path.close();
	//// cout << path_count << " path in double direction method" << endl;
	return ;

}

paths find_all_paths_use_parents(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  src, NODE_TYPE  dst, NODE_TYPE  distance, map<NODE_TYPE, vector<NODE_TYPE > >* parents)
{
	paths result;
	//output all the path
	stack<pair<NODE_TYPE, NODE_TYPE > > path;
	NODE_TYPE cur_distance = distance;


	for (cur_distance = 1; cur_distance <= distance; cur_distance++)
	{
		paths temp = find_all_paths_use_parents_single_distance(adjacency_list, node_num, src, dst, cur_distance, parents);
		result.add_paths(temp);
	}
	return result;
}
 

paths find_all_meet_paths(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num,NODE_TYPE src, NODE_TYPE dst, meet_nodes nodes_nohot, map<NODE_TYPE, vector<NODE_TYPE > >* parents_left, map<NODE_TYPE, vector<NODE_TYPE > >*parents_right, NODE_TYPE k)
{
	paths result;
	//output all the path
	stack<pair<NODE_TYPE, NODE_TYPE > > path;
	NODE_TYPE meet_node, left_distance, right_distandce;
	for (auto iter = nodes_nohot.get_nodes_pairs().begin(); iter != nodes_nohot.get_nodes_pairs().end(); iter++)
	{
		meet_node = iter->node;
		left_distance = iter->left_distance;
		right_distandce = iter->right_distance;
		if (left_distance + right_distandce > k) { continue; }
		//// cout << meet_node << " : " << left_distance << " : " << right_distandce << " in meet path" << endl;
		
		paths left_part = find_all_paths_use_parents_single_distance(adjacency_list, node_num, src, meet_node, left_distance,parents_left);
		paths right_part = find_all_paths_use_parents_single_distance(adjacency_list_reverse, node_num, dst, meet_node, right_distandce, parents_right);
		right_part.reverse();
		paths temp = left_part.join(right_part,k);
		result.add_paths(temp);
	}
	return result;
}




paths find_all_paths_use_parents_single_distance(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  src, NODE_TYPE  dst, NODE_TYPE  distance, map<NODE_TYPE, vector<NODE_TYPE > >* parents)
{
	paths result;
	if (distance == 0)
	{
		vector<NODE_TYPE> src_path;
		src_path.push_back(src);
		result.push_back(src_path);
		return result;
	}
	//output all the path
	stack<pair<NODE_TYPE, NODE_TYPE > > path;
	NODE_TYPE cur_distance = distance;


	auto iter4 = parents[dst].find(cur_distance);
	if (iter4 == parents[dst].end()) {//if no vector
		return result;
	}
	for (auto iter = iter4->second.begin(); iter != iter4->second.end(); iter++)
	{
		path.push(pair<NODE_TYPE, NODE_TYPE >(*iter, cur_distance));
	}

	vector<NODE_TYPE> output_path;
	output_path.push_back(dst);
	while (!path.empty())
	{
		NODE_TYPE  temp_node = path.top().first;
		NODE_TYPE  distance = path.top().second;
		if (temp_node == src && distance == 1)
		{
			//output_path
			if (!circle_detection(output_path))
			{
				vector<NODE_TYPE> temp = output_path;
				temp.push_back(src);
				reverse(temp.begin(), temp.end());
				result.push_back(temp);
			}
			path.pop();
			continue;
		}
		else if (temp_node == output_path.back())
		{
			output_path.pop_back();
			path.pop();
			continue;
		}
		else {
			output_path.push_back(temp_node);
		}

		distance--;
		auto iter = parents[temp_node].find(distance);
		auto iter2 = iter->second;
		for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
		{
			path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
		}
	}
	return result;
}



paths find_all_paths_with_hotpoints(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points, path_index index)
{
	map<NODE_TYPE, vector<NODE_TYPE > >* parents_left = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];
	map<NODE_TYPE, vector<NODE_TYPE > >* parents_right = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];
	map<NODE_TYPE, vector<NODE_TYPE> > left_hot;// distance,vector(node)
	map<NODE_TYPE, vector<NODE_TYPE> > right_hot;
	paths paths_without_hot, hot_paths;
	//NODE_TYPE cur_distance = 1;
	NODE_TYPE left_distance = 0;
	NODE_TYPE right_distance = 0;
	set<NODE_TYPE> left_propagate;
	set<NODE_TYPE> right_propagate;
	left_propagate.insert(query_node1);
	right_propagate.insert(query_node2);
	while (left_distance < k)
	{
		if (hot_points.find(query_node1) != hot_points.end())//judge whether query node 1 is a hot point, if yes, we stop left traversal
		{
			vector<NODE_TYPE> a;
			a.push_back(query_node1);
			left_hot.insert(pair<NODE_TYPE, vector<NODE_TYPE > >(0, a));
			break;
		}
		left_distance += 1;
		set<NODE_TYPE > temp_propagate;
		for (auto iter = left_propagate.begin(); iter != left_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes
				if (parents_left[*iter2].find(left_distance) == parents_left[*iter2].end()) {
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents_left[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(left_distance, a));
				}
				else {
					auto iter3 = parents_left[*iter2].find(left_distance);
					iter3->second.push_back(cur_node);//insert node distance value
				}


				//stop condition
				if (hot_points.find(*iter2) != hot_points.end())//*iter2 is a hot points
				{
					// add this hot points into left hot vector
					if (left_hot.find(*iter2) == left_hot.end())//this node is the first to insert
					{
						vector<NODE_TYPE> a;
						a.push_back(*iter2);
						left_hot.insert(pair<NODE_TYPE, vector<NODE_TYPE > >(left_distance, a));
					}
					else {
						auto iter3 = left_hot.find(left_distance);
						iter3->second.push_back(*iter2);
					}
					continue;
				}
				//end of stop condition

				temp_propagate.insert(*iter2);//insert all the new propagate nodes
			}
		}
		left_propagate = temp_propagate;
	}
	while (right_distance < k)
	{

		if (hot_points.find(query_node2) != hot_points.end())
		{
			vector<NODE_TYPE> a;
			a.push_back(query_node2);
			right_hot.insert(pair<NODE_TYPE, vector<NODE_TYPE > >(0, a));
			break;
		}
		right_distance += 1;
		set<NODE_TYPE > temp_propagate;
		for (auto iter = right_propagate.begin(); iter != right_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			for (auto iter2 = adjacency_list_reverse[cur_node].begin(); iter2 != adjacency_list_reverse[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes
				if (parents_right[*iter2].find(right_distance) == parents_right[*iter2].end()) {
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents_right[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(right_distance, a));
				}
				else {
					auto iter3 = parents_right[*iter2].find(right_distance);
					iter3->second.push_back(cur_node);//insert node distance value
				}
				//stop condition
				if (hot_points.find(*iter2) != hot_points.end())//*iter2 is a hot points
				{
					// add this hot points into left hot vector
					if (right_hot.find(*iter2) == right_hot.end())//this node is the first to insert
					{
						vector<NODE_TYPE> a;
						a.push_back(*iter2);
						right_hot.insert(pair<NODE_TYPE, vector<NODE_TYPE > >(right_distance, a));
					}
					else {
						auto iter3 = right_hot.find(right_distance);
						iter3->second.push_back(*iter2);
					}
					continue;
				}
				//end of stop condition
				temp_propagate.insert(*iter2);//insert all the new propagate nodes
			}
		}
		right_propagate = temp_propagate;
	}
	paths_without_hot = find_all_paths_use_parents(adjacency_list, node_num, query_node1, query_node2, k, parents_left);

	// cout << left_hot.size() << ":" << right_hot.size() << endl;
	//to find hot_paths, judge whether there is a hot path first
	NODE_TYPE left_node, right_node, temp_left_distance, temp_right_distance;
	//NODE_TYPE left_min_dis, right_min_dis;
	for (auto iter1 = left_hot.begin(); iter1 != left_hot.end(); iter1++)
	{
		for (auto iter2 = right_hot.begin(); iter2 != right_hot.end(); iter2++)
		{
			paths left_hot_pahts_single_node, right_hot_pahts_single_node;
			paths temp_left;
			paths temp_right;
			for (auto iter_left_node = iter1->second.begin(); iter_left_node != iter1->second.end(); iter_left_node++)
			{
				for (auto iter_right_node = iter2->second.begin(); iter_right_node != iter2->second.end(); iter_right_node++)
				{
					left_node = *iter_left_node;
					right_node = *iter_right_node;
					temp_left_distance = iter1->first;
					temp_right_distance = iter2->first;
					// cout << left_node << " : " << right_node << " : " << temp_left_distance << " : " << temp_right_distance << " : " << endl;
					temp_left = find_all_paths_use_parents_single_distance(adjacency_list, node_num, query_node1, left_node, temp_left_distance, parents_left);
					temp_right = find_all_paths_use_parents_single_distance(adjacency_list_reverse, node_num, query_node2, right_node, temp_right_distance, parents_right);
					temp_right.reverse();
					left_hot_pahts_single_node.add_paths(temp_left);
					right_hot_pahts_single_node.add_paths(temp_right);
					paths paths_bt_hot_points;
					paths_bt_hot_points = index.find_paths_between_two_hot_nodes(left_node, right_node, (k - temp_left_distance - temp_right_distance));
					paths temp_result;
					// cout << " left temp:";
					temp_left.output();
					// cout << " right temp:";
					temp_right.output();
					// cout << " paths bt hot points :";
					paths_bt_hot_points.output();
					temp_result = temp_left.join(paths_bt_hot_points, k);
					temp_result = temp_result.join(temp_right, k);
					hot_paths.add_paths(temp_result);
				}
			}

		}
	}
	paths result;
	hot_paths.drop_path_length_less_than_k(1);// hot points path may contains path with length equal to 1
	result.add_paths(paths_without_hot);
	result.add_paths(hot_paths);
	hot_paths.write_to_file("hint.txt");
	return result;
}



paths find_all_paths_with_hotpoints_step_by_step(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points, path_index &index)
{
	map<NODE_TYPE, vector<NODE_TYPE > >* parents_left = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];
	map<NODE_TYPE, vector<NODE_TYPE > >* parents_right = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];
	map<NODE_TYPE, vector<NODE_TYPE> > left_hot;// distance,vector(node)
	map<NODE_TYPE, vector<NODE_TYPE> > right_hot;
	paths paths_without_hot, hot_paths;
	NODE_TYPE left_distance = 0;
	NODE_TYPE right_distance = 0;
	set<NODE_TYPE> left_propagate;
	set<NODE_TYPE> right_propagate;
	left_propagate.insert(query_node1);
	right_propagate.insert(query_node2);
	vector<NODE_TYPE > temp_n1, temp_n2;
	temp_n1.push_back(query_node1);
	temp_n2.push_back(query_node2);
	parents_left[query_node1].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(0, temp_n1));
	parents_left[query_node2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(0, temp_n2));

	bool left_go_on = true;
	bool right_go_on = true;

	NODE_TYPE left_hot_point_distance = 0;
	NODE_TYPE right_hot_point_distance = 0;
	bool left_hot_distance_updated, rihgt_hot_distance_updated;
	left_hot_distance_updated = false;
	rihgt_hot_distance_updated = false;
	// some structures to store unhot paths
	meet_nodes nodes_unhot;
	while (left_go_on || right_go_on)
	{
		if (left_distance < (k-right_hot_point_distance) && left_go_on)
		{
			//// cout << " left distance: " << left_distance << endl;
			if (hot_points.find(query_node1) != hot_points.end())//judge whether query node 1 is a hot point, if yes, we stop left travesal
			{
				if (!left_hot_distance_updated)
				{
					left_hot_distance_updated = true;
					left_hot_point_distance = left_distance;
				}
				vector<NODE_TYPE> a;
				a.push_back(query_node1);
				left_hot.insert(pair<NODE_TYPE, vector<NODE_TYPE > >(0, a));
				left_go_on = false;
			}
			else
			{
				left_distance += 1;
				set<NODE_TYPE > temp_propagate;
				for (auto iter = left_propagate.begin(); iter != left_propagate.end(); iter++)//traverse all the cur_propagate set nodes
				{//insert <node,distance> pair first
					NODE_TYPE  cur_node = *iter;
					for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
					{//traverse all the adjacent nodes
						if (parents_right[*iter2].find(left_distance - 1) != parents_right[*iter2].end()) //there is a right meet node
						{
							auto right_meet = parents_right[*iter2].find(left_distance - 1);
							meet_node a(*iter2, left_distance, left_distance - 1);
							nodes_unhot.push_back(a);
						}
						if (parents_left[*iter2].find(left_distance) == parents_left[*iter2].end()) {
							vector<NODE_TYPE > a;
							a.push_back(cur_node);
							parents_left[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(left_distance, a));
						}
						else {
							auto iter3 = parents_left[*iter2].find(left_distance);
							iter3->second.push_back(cur_node);//insert node distance value
						}


						//stop condition
						if (hot_points.find(*iter2) != hot_points.end())//*iter2 is a hot points
						{
							if (!left_hot_distance_updated)
							{
								left_hot_distance_updated = true;
								left_hot_point_distance = left_distance;
							}
							// add this hot points into left hot vector
							if (left_hot.find(*iter2) == left_hot.end())//this node is the first to insert
							{
								vector<NODE_TYPE> a;
								a.push_back(*iter2);
								left_hot.insert(pair<NODE_TYPE, vector<NODE_TYPE > >(left_distance, a));
							}
							else {
								auto iter3 = left_hot.find(left_distance);
								iter3->second.push_back(*iter2);
							}
							continue;
						}
						//end of stop condition

						temp_propagate.insert(*iter2);//insert all the new propagate nodes
					}
				}
				left_propagate = temp_propagate;
			}
		}
		else { left_go_on = false; }

		if (right_distance < (k - left_hot_point_distance) && right_go_on)
		{
			//// cout << " right distance: " << right_distance << endl;
			if (hot_points.find(query_node2) != hot_points.end())
			{
				if (!rihgt_hot_distance_updated)
				{
					rihgt_hot_distance_updated = true;
					right_hot_point_distance = right_distance;
				}
				vector<NODE_TYPE> a;
				a.push_back(query_node2);
				right_hot.insert(pair<NODE_TYPE, vector<NODE_TYPE > >(0, a));
				right_go_on = false;
			}
			else
			{
				right_distance += 1;
				set<NODE_TYPE > temp_propagate;
				for (auto iter = right_propagate.begin(); iter != right_propagate.end(); iter++)//traverse all the cur_propagate set nodes
				{//insert <node,distance> pair first
					NODE_TYPE  cur_node = *iter;
					for (auto iter2 = adjacency_list_reverse[cur_node].begin(); iter2 != adjacency_list_reverse[cur_node].end(); iter2++)
					{//traverse all the adjacent nodes
						if (parents_left[*iter2].find(right_distance) != parents_left[*iter2].end()) //there is a left meet node
						{
							auto left_meet = parents_left[*iter2].find(right_distance);
							meet_node a(*iter2, right_distance, right_distance);
							nodes_unhot.push_back(a);
						}


						if (parents_right[*iter2].find(right_distance) == parents_right[*iter2].end()) {
							vector<NODE_TYPE > a;
							a.push_back(cur_node);
							parents_right[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(right_distance, a));
						}
						else {
							auto iter3 = parents_right[*iter2].find(right_distance);
							iter3->second.push_back(cur_node);//insert node distance value
						}
						//stop condition
						if (hot_points.find(*iter2) != hot_points.end())//*iter2 is a hot points
						{
							if (!rihgt_hot_distance_updated)
							{
								rihgt_hot_distance_updated = true;
								right_hot_point_distance = right_distance;
							}
							// add this hot points into left hot vector
							if (right_hot.find(*iter2) == right_hot.end())//this node is the first to insert
							{
								vector<NODE_TYPE> a;
								a.push_back(*iter2);
								right_hot.insert(pair<NODE_TYPE, vector<NODE_TYPE > >(right_distance, a));
							}
							else {
								auto iter3 = right_hot.find(right_distance);
								iter3->second.push_back(*iter2);
							}
							continue;
						}
						//end of stop condition
						temp_propagate.insert(*iter2);//insert all the new propagate nodes
					}
				}
				right_propagate = temp_propagate;
			}
		}
		else { right_go_on = false; }
	}
	// cout << " size of meetnodes : " << nodes_unhot.get_nodes_pairs().size() << endl;
	paths_without_hot = find_all_meet_paths(adjacency_list,adjacency_list_reverse,node_num,query_node1,query_node2,nodes_unhot,parents_left,parents_right,k);
	//paths_without_hot.output();
	// cout << left_hot.size() << ":" << right_hot.size() << endl;
	//to find hot_paths, judge whether there is a hot path first
	NODE_TYPE left_node, right_node, temp_left_distance, temp_right_distance;
	paths index_new_paths;
	for (auto iter1 = left_hot.begin(); iter1 != left_hot.end(); iter1++)
	{
		//// cout << " iter 1" << endl;
		for (auto iter2 = right_hot.begin(); iter2 != right_hot.end(); iter2++)
		{
			//// cout << "iter 2" << endl;
			paths left_hot_pahts_single_node, right_hot_pahts_single_node;
			paths temp_left;
			paths temp_right;
			for (auto iter_left_node = iter1->second.begin(); iter_left_node != iter1->second.end(); iter_left_node++)
			{
				for (auto iter_right_node = iter2->second.begin(); iter_right_node != iter2->second.end(); iter_right_node++)
				{
					left_node = *iter_left_node;
					right_node = *iter_right_node;
					temp_left_distance = iter1->first;
					temp_right_distance = iter2->first;
					//// cout << left_node << " : " << right_node << " : " << temp_left_distance << " : " << temp_right_distance << " : " << endl;
					temp_left = find_all_paths_use_parents_single_distance(adjacency_list, node_num, query_node1, left_node, temp_left_distance, parents_left);
					temp_right = find_all_paths_use_parents_single_distance(adjacency_list_reverse, node_num, query_node2, right_node, temp_right_distance, parents_right);
					temp_right.reverse();
					left_hot_pahts_single_node.add_paths(temp_left);
					right_hot_pahts_single_node.add_paths(temp_right);
					paths paths_bt_hot_points;
					cur_path c_path;
					index.find_paths_between_two_hot_nodes_index_without_cross_other_hotpoints(paths_bt_hot_points, c_path,left_node, right_node,(k - temp_left_distance - temp_right_distance),0);

					paths temp_result;

					//// cout << " left temp:";
					//temp_left.output();
					//// cout << " right temp:";
					//temp_right.output();
					//// cout << " paths bt hot points :";
					//paths_bt_hot_points.output();
					temp_result = temp_left.join(paths_bt_hot_points, k);

					temp_result = temp_result.join(temp_right, k);

					hot_paths.add_paths(temp_result);
					
					//update the index use known path, this dynamic update is not correct, it miss some path. So we comment it.s
					paths index_new;
					vector<vector<NODE_TYPE> > edge;
					vector<NODE_TYPE> query_edge;
					query_edge.push_back(query_node2);
					query_edge.push_back(query_node1);// the new updated edge is from query_node2 to query_node1(the direction)
					edge.push_back(query_edge);
					paths updated_edge(edge);
					index_new = temp_right.join(updated_edge,k);
					index_new = index_new.join(temp_left,k);
					index_new_paths.add_paths(index_new);
					
					//end of update
				}
			}

		}
	}
	index.push_back(index_new_paths);
	paths result;
	hot_paths.drop_path_length_less_than_k(1);// hot points path may contains path with length equal to 1
	result.add_paths(paths_without_hot);
	result.add_paths(hot_paths);
	hot_paths.write_to_file("hint_step_by_step.txt");
	return result;

}


paths find_all_paths_with_hotpoints_dfs(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points, path_index &index)
{
	int min_stop_distance = int(k);
	int right_stop_distance = int(k);

	// some structures to store unhot paths
	set<NODE_TYPE> stop_points(hot_points);
	stop_points.insert(query_node2);
	map<NODE_TYPE, paths> left_dfs;
	map<NODE_TYPE, paths> right_dfs;
	vector<NODE_TYPE> c_path;
	set<NODE_TYPE> c_path_set;
	dfs(adjacency_list_double, adjacency_list, left_dfs, query_node1, c_path, stop_points, k, 0, min_stop_distance, c_path_set);

	paths paths_without_hot;
	auto iter_left = left_dfs.find(query_node2);
	if (iter_left != left_dfs.end())
	{
		paths_without_hot.add_paths(iter_left->second);
		//paths_without_hot.output();
	}
	//paths_without_hot.output();
	set<NODE_TYPE> right_stop_points(hot_points);
	right_stop_points.insert(query_node1);
	dfs(adjacency_list_double, adjacency_list_reverse, right_dfs, query_node2, c_path, right_stop_points, k-min_stop_distance, 0, right_stop_distance, c_path_set);

	//to find hot_paths, judge whether there is a hot path first
	paths hot_paths;
	NODE_TYPE left_node, right_node, temp_left_distance, temp_right_distance;
	paths index_new_paths;
	//// cout << " start dfs find hot points " << endl;
	for (auto iter1 = left_dfs.begin(); iter1 != left_dfs.end(); iter1++)
	{
		if (iter1->first == query_node2)// we should remove query_node2 from left_dfs hot_points to find paths_bt_hot_points
		{
			continue;
		}
		for (auto iter2 = right_dfs.begin(); iter2 != right_dfs.end(); iter2++)
		{
			paths temp_result;
			paths paths_bt_hot_points;

			//// cout << "iter 2" << endl;
			//paths left_hot_pahts_single_node, right_hot_pahts_single_node;
			paths temp_left;
			paths temp_right;
			left_node = iter1->first;
			right_node = iter2->first;
			//// cout << left_node << " left: " << right_node << " right: " << endl;// << temp_left_distance << " : " << temp_right_distance << " : " << endl;
			temp_left = iter1->second;
			temp_right = iter2->second;
			temp_right.reverse();
			if (left_node == right_node)
			{
				temp_result = temp_left.join(paths_bt_hot_points, k);
				temp_result = temp_result.join(temp_right, k);


				hot_paths.add_paths(temp_result);
				continue;
			}
			//left_hot_pahts_single_node.add_paths(temp_left);
			//right_hot_pahts_single_node.add_paths(temp_right);
			temp_left_distance = temp_left.get_min_distance();
			temp_right_distance = temp_right.get_min_distance();
			cur_path c_path;
			index.find_paths_between_two_hot_nodes_index_without_cross_other_hotpoints(paths_bt_hot_points, c_path, left_node, right_node, (k - temp_left_distance - temp_right_distance), 0);
			//// cout << " paths bt hot points " << endl;
			//paths_bt_hot_points.output();
			

			//// cout << " left temp:";
			//temp_left.output();
			//// cout << " right temp:";
			//temp_right.output();
			//// cout << " paths bt hot points :";
			//paths_bt_hot_points.output();
			temp_result = temp_left.join(paths_bt_hot_points, k);
			temp_result = temp_result.join(temp_right, k);


			hot_paths.add_paths(temp_result);

			//update the index use known path, this dynamic update is not correct, it miss some path. So we comment it.s
			paths index_new;
			vector<vector<NODE_TYPE> > edge;
			vector<NODE_TYPE> query_edge;
			query_edge.push_back(query_node2);
			query_edge.push_back(query_node1);// the new updated edge is from query_node2 to query_node1(the direction)
			edge.push_back(query_edge);
			paths updated_edge(edge);
			index_new = temp_right.join(updated_edge, k);
			index_new = index_new.join(temp_left, k);
			index_new_paths.add_paths(index_new);
			
			//end of update	

		}
	}
	bool updateIndex = true;
	//if the new edge is already in the adjacency_list, then we do not need to update index
	for (auto iter = adjacency_list[query_node2].begin(); iter != adjacency_list[query_node2].end(); iter++)
	{
		if (*iter == query_node1) {
			updateIndex = false;
			break;
		}
	}
	if (updateIndex) {
		index_new_paths.drop_path_length_less_than_k(1);
		index.push_back(index_new_paths);//already updated in hp index algorithm
		
		//update index
		//index.push_back(index_new_paths);
		//adjacency_list[query_node2].push_back(query_node1);//if the query is from node1 to node2, then the new edge is node2->node1 to from a circle
		//adjacency_list_reverse[query_node1].push_back(query_node2);
		//since adjacency_list and index.push back will be updated in the other algorithm
		//end of update
	}
	paths result;
	hot_paths.drop_path_length_less_than_k(1);// hot points path may contains path with length equal to 1
	result.add_paths(paths_without_hot);
	result.add_paths(hot_paths);
	result.drop_path_with_repeat_node();
	result.drop_path_not_start_from_nodeS(query_node1);
	//result.sort_by_string_order();
	//hot_paths.write_to_file("hint_step_by_step.txt");
	return result;

}

paths find_all_paths_with_hotpoints_dfs_get_update_time_new_algor(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, unordered_map<NODE_TYPE, set<NODE_TYPE>>& induced_subgraph, unordered_map<NODE_TYPE, set<NODE_TYPE>>& induced_subgraph_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points, path_index &index, double& update_time)
{
	update_time = 0;
	clock_t start, end;
	int min_stop_distance = int(k);
	int right_stop_distance = int(k);

	// some structures to store unhot paths
	set<NODE_TYPE> stop_points(hot_points);
	stop_points.insert(query_node2);
	map<NODE_TYPE, paths> left_dfs;
	map<NODE_TYPE, paths> right_dfs;
	vector<NODE_TYPE> c_path;
	set<NODE_TYPE> c_path_set;

	pruned_subgraph_unordered_map_double_direction temp_pruned;
	temp_pruned.meet_nodes = stop_points;
	left_dfs = temp_pruned.dfs_without_recursion_all_meetpoints(induced_subgraph, query_node1, NULL_NODE, k, node_num);

	//left_dfs = ;
	//dfs_ex_node2(adjacency_list_double, adjacency_list, left_dfs, query_node1, c_path, stop_points, k, 0, min_stop_distance, c_path_set, query_node2);

	paths paths_without_hot;
	auto iter_left = left_dfs.find(query_node2);
	if (iter_left != left_dfs.end())
	{
		if(iter_left->second.path.size()-1 < min_stop_distance)
		{
			min_stop_distance = iter_left->second.path.size() - 1;
		}
		paths_without_hot.add_paths(iter_left->second);
		//paths_without_hot.output();
	}
	if(paths_without_hot.path.size() == 0 )
	{
		return paths_without_hot;
	}
	//paths_without_hot.output();
	set<NODE_TYPE> right_stop_points(hot_points);
	right_stop_points.insert(query_node1);
	temp_pruned.meet_nodes = right_stop_points;
	right_dfs = temp_pruned.dfs_without_recursion_all_meetpoints(induced_subgraph_reverse, query_node2, NULL_NODE, k-min_stop_distance, node_num);
	//dfs_ex_node2(adjacency_list_double, adjacency_list_reverse, right_dfs, query_node2, c_path, right_stop_points, k - min_stop_distance, 0, right_stop_distance, c_path_set, query_node1);

	//to find hot_paths, judge whether there is a hot path first
	paths hot_paths;
	NODE_TYPE left_node, right_node, temp_left_distance, temp_right_distance;
	paths index_new_paths;
	//// cout << " start dfs find hot points " << endl;
	for (auto iter1 = left_dfs.begin(); iter1 != left_dfs.end(); iter1++)
	{
		if (iter1->first == query_node2)// we should remove query_node2 from left_dfs hot_points to find paths_bt_hot_points
		{
			continue;
		}
		for (auto iter2 = right_dfs.begin(); iter2 != right_dfs.end(); iter2++)
		{
			paths temp_result;
			paths paths_bt_hot_points;

			//// cout << "iter 2" << endl;
			//paths left_hot_pahts_single_node, right_hot_pahts_single_node;
			paths temp_left;
			paths temp_right;
			left_node = iter1->first;
			right_node = iter2->first;
			//// cout << left_node << " left: " << right_node << " right: " << endl;// << temp_left_distance << " : " << temp_right_distance << " : " << endl;
			temp_left = iter1->second;
			temp_right = iter2->second;
			temp_right.reverse();
			if (left_node == right_node)
			{
				temp_result = temp_left.join(paths_bt_hot_points, k);
				temp_result = temp_result.join(temp_right, k);


				hot_paths.add_paths(temp_result);
				continue;
			}
			//left_hot_pahts_single_node.add_paths(temp_left);
			//right_hot_pahts_single_node.add_paths(temp_right);
			temp_left_distance = temp_left.get_min_distance();
			temp_right_distance = temp_right.get_min_distance();
			cur_path c_path;
			index.find_paths_between_two_hot_nodes_index_without_cross_other_hotpoints(paths_bt_hot_points, c_path, left_node, right_node, (k - temp_left_distance - temp_right_distance), 0);
			//// cout << " paths bt hot points " << endl;
			//paths_bt_hot_points.output();


			//// cout << " left temp:";
			//temp_left.output();
			//// cout << " right temp:";
			//temp_right.output();
			//// cout << " paths bt hot points :";
			//paths_bt_hot_points.output();
			temp_result = temp_left.join(paths_bt_hot_points, k);
			temp_result = temp_result.join(temp_right, k);


			hot_paths.add_paths(temp_result);
			start = clock();
			//update the index use known path, this dynamic update is not correct, it miss some path. So we comment it.s
			paths index_new;
			vector<vector<NODE_TYPE> > edge;
			vector<NODE_TYPE> query_edge;
			query_edge.push_back(query_node2);
			query_edge.push_back(query_node1);// the new updated edge is from query_node2 to query_node1(the direction)
			edge.push_back(query_edge);
			paths updated_edge(edge);
			index_new = temp_right.join(updated_edge, k);
			index_new = index_new.join(temp_left, k);
			index_new_paths.add_paths(index_new);
			end = clock();
			update_time += double(end - start) / CLOCKS_PER_SEC;
			//end of update	

		}
	}
	bool updateIndex = false;
	//if the new edge is already in the adjacency_list, then we do not need to update index
	for (auto iter = adjacency_list[query_node2].begin(); iter != adjacency_list[query_node2].end(); iter++)
	{
		if (*iter == query_node1) {
			updateIndex = false;
			break;
		}
	}
	if (updateIndex) {
		index_new_paths.drop_path_length_less_than_k(1);
		index.push_back(index_new_paths);//already updated in hp index algorithm

										 //update index
										 //index.push_back(index_new_paths);
										 //adjacency_list[query_node2].push_back(query_node1);//if the query is from node1 to node2, then the new edge is node2->node1 to from a circle
										 //adjacency_list_reverse[query_node1].push_back(query_node2);
										 //since adjacency_list and index.push back will be updated in the other algorithm
										 //end of update
	}
	paths result;
	hot_paths.drop_path_length_less_than_k(1);// hot points path may contains path with length equal to 1
	result.add_paths(paths_without_hot);
	result.add_paths(hot_paths);
	//result.sort_by_string_order();
	//hot_paths.write_to_file("hint_step_by_step.txt");
	return result;

}



paths find_all_paths_with_hotpoints_dfs_get_update_time(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points, path_index &index, double& update_time)
{
	update_time = 0;
	clock_t start, end;
	int min_stop_distance = int(k);
	int right_stop_distance = int(k);

	// some structures to store unhot paths
	set<NODE_TYPE> stop_points(hot_points);
	stop_points.insert(query_node2);
	map<NODE_TYPE, paths> left_dfs;
	map<NODE_TYPE, paths> right_dfs;
	vector<NODE_TYPE> c_path;
	set<NODE_TYPE> c_path_set;
	dfs_ex_node2(adjacency_list_double, adjacency_list, left_dfs, query_node1, c_path, stop_points, k, 0, min_stop_distance, c_path_set,query_node2);

	paths paths_without_hot;
	auto iter_left = left_dfs.find(query_node2);
	if (iter_left != left_dfs.end())
	{
		paths_without_hot.add_paths(iter_left->second);
		//paths_without_hot.output();
	}
	//paths_without_hot.output();
	set<NODE_TYPE> right_stop_points(hot_points);
	right_stop_points.insert(query_node1);
	dfs_ex_node2(adjacency_list_double, adjacency_list_reverse, right_dfs, query_node2, c_path, right_stop_points, k - min_stop_distance, 0, right_stop_distance, c_path_set, query_node1);

	//to find hot_paths, judge whether there is a hot path first
	paths hot_paths;
	NODE_TYPE left_node, right_node, temp_left_distance, temp_right_distance;
	paths index_new_paths;
	//// cout << " start dfs find hot points " << endl;
	for (auto iter1 = left_dfs.begin(); iter1 != left_dfs.end(); iter1++)
	{
		if (iter1->first == query_node2)// we should remove query_node2 from left_dfs hot_points to find paths_bt_hot_points
		{
			continue;
		}
		for (auto iter2 = right_dfs.begin(); iter2 != right_dfs.end(); iter2++)
		{
			paths temp_result;
			paths paths_bt_hot_points;

			//// cout << "iter 2" << endl;
			//paths left_hot_pahts_single_node, right_hot_pahts_single_node;
			paths temp_left;
			paths temp_right;
			left_node = iter1->first;
			right_node = iter2->first;
			//// cout << left_node << " left: " << right_node << " right: " << endl;// << temp_left_distance << " : " << temp_right_distance << " : " << endl;
			temp_left = iter1->second;
			temp_right = iter2->second;
			temp_right.reverse();
			if (left_node == right_node)
			{
				temp_result = temp_left.join(paths_bt_hot_points, k);
				temp_result = temp_result.join(temp_right, k);


				hot_paths.add_paths(temp_result);
				continue;
			}
			//left_hot_pahts_single_node.add_paths(temp_left);
			//right_hot_pahts_single_node.add_paths(temp_right);
			temp_left_distance = temp_left.get_min_distance();
			temp_right_distance = temp_right.get_min_distance();
			cur_path c_path;
			index.find_paths_between_two_hot_nodes_index_without_cross_other_hotpoints(paths_bt_hot_points, c_path, left_node, right_node, (k - temp_left_distance - temp_right_distance), 0);
			//// cout << " paths bt hot points " << endl;
			//paths_bt_hot_points.output();


			//// cout << " left temp:";
			//temp_left.output();
			//// cout << " right temp:";
			//temp_right.output();
			//// cout << " paths bt hot points :";
			//paths_bt_hot_points.output();
			temp_result = temp_left.join(paths_bt_hot_points, k);
			temp_result = temp_result.join(temp_right, k);


			hot_paths.add_paths(temp_result);
			start = clock();
			//update the index use known path, this dynamic update is not correct, it miss some path. So we comment it.s
			paths index_new;
			vector<vector<NODE_TYPE> > edge;
			vector<NODE_TYPE> query_edge;
			query_edge.push_back(query_node2);
			query_edge.push_back(query_node1);// the new updated edge is from query_node2 to query_node1(the direction)
			edge.push_back(query_edge);
			paths updated_edge(edge);
			index_new = temp_right.join(updated_edge, k);
			index_new = index_new.join(temp_left, k);
			index_new_paths.add_paths(index_new);
			end = clock();
			update_time += double(end - start) / CLOCKS_PER_SEC;
			//end of update	

		}
	}
	bool updateIndex = true;
	//if the new edge is already in the adjacency_list, then we do not need to update index
	for (auto iter = adjacency_list[query_node2].begin(); iter != adjacency_list[query_node2].end(); iter++)
	{
		if (*iter == query_node1) {
			updateIndex = false;
			break;
		}
	}
	if (updateIndex) {
		index_new_paths.drop_path_length_less_than_k(1);
		index.push_back(index_new_paths);//already updated in hp index algorithm

										 //update index
										 //index.push_back(index_new_paths);
										 //adjacency_list[query_node2].push_back(query_node1);//if the query is from node1 to node2, then the new edge is node2->node1 to from a circle
										 //adjacency_list_reverse[query_node1].push_back(query_node2);
										 //since adjacency_list and index.push back will be updated in the other algorithm
										 //end of update
	}
	paths result;
	hot_paths.drop_path_length_less_than_k(1);// hot points path may contains path with length equal to 1
	result.add_paths(paths_without_hot);
	result.add_paths(hot_paths);
	//result.sort_by_string_order();
	//hot_paths.write_to_file("hint_step_by_step.txt");
	return result;

}



paths single_direction_baseline_stop_at_hotpoints_2(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points)
{//original construct index method
	set<NODE_TYPE> visited_node;
	visited_node.insert(query_node2);
	visited_node.insert(query_node1);
	paths result_path;
	//ofstream outPath;
	//outPath.open("out.txt");
	NODE_TYPE path_count = 0;
	map<NODE_TYPE, vector<NODE_TYPE > >* parents = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];
	set<NODE_TYPE > cur_propagate;
	NODE_TYPE  cur_distance = 1;
	cur_propagate.insert(query_node1);
	while (cur_distance <= k)
	{
		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes
				if (visited_node.find(*iter2) != visited_node.end())
				{
					continue;
				}
				if (*iter2 != query_node2 && (hot_points.find(*iter2) != hot_points.end()))
				{
					continue;
				}
				if (parents[*iter2].find(cur_distance) == parents[*iter2].end()) {
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
				}
				if (*iter2 == query_node2)
				{
				}
				else
				{

					temp_propagate.insert(*iter2);//insert all the new propagate nodes
					visited_node.insert(*iter2);
				}
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}
	//output all the path
	result_path = find_all_paths_use_parents(adjacency_list, node_num, query_node1, query_node2, k, parents);
	delete[] parents;
	return result_path;
}

// to code here 
paths single_direction_baseline_stop_at_hotpoints_dfs(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points)
{//original construct index method
	paths result_path;
	stack<NODE_TYPE> c_path;
	vector<NODE_TYPE> c_path_v;
	set<NODE_TYPE> c_path_set;
	int cur_distance = 1;
	vector< vector<NODE_TYPE>::iterator > cur_iters;
	NODE_TYPE src = query_node1;
	c_path.push(src);
	c_path_v.push_back(src);
	vector<int> min_order;
	for (int i = 0; i <= k + 1; i++)
	{
		cur_iters.push_back(adjacency_list[src].begin());
	}
	//set<NODE_TYPE> new_c_path_set(c_path_set);
	c_path_set.insert(src);
	bool dfs_finish = false;
	NODE_TYPE cur_node;
	while (!dfs_finish)
	{// we do not use adjacency_list_double
		NODE_TYPE cur_node, next_node;
		if (cur_distance <= 0)
		{
			break;
		}
		auto iter = cur_iters[cur_distance];
		if (cur_distance == 1)
		{
			if (cur_iters[cur_distance] == adjacency_list[src].end())
			{
				cur_distance--;
				cur_iters[cur_distance]++;
				cur_node = c_path.top();
				c_path.pop();
				c_path_v.pop_back();
				auto iter2 = c_path_set.find(cur_node);
				if (iter2 != c_path_set.end())
				{
					c_path_set.erase(iter2);
				}
				continue;
			}
		}
		else if (cur_iters[cur_distance] == adjacency_list[*cur_iters[cur_distance - 1]].end() || cur_distance > k)
		{
			// go back
			cur_distance--;
			cur_iters[cur_distance]++;
			cur_node = c_path.top();
			c_path.pop();
			c_path_v.pop_back();
			auto iter2 = c_path_set.find(cur_node);
			if (iter2 != c_path_set.end())
			{
				c_path_set.erase(iter2);
			}
			continue;
		}
		cur_node = *iter;
		//normal logic
		bool force_continue = false;
		while (c_path_set.find(*cur_iters[cur_distance]) != c_path_set.end())//*cur_iters[cur_distance] < src || 
		{
			cur_iters[cur_distance]++;
			if (cur_distance == 1)
			{
				if (cur_iters[cur_distance] == adjacency_list[src].end())
				{
					force_continue = true;
					break;
				}
			}
			else if (cur_iters[cur_distance] == adjacency_list[*cur_iters[cur_distance - 1]].end())
			{
				force_continue = true;
				break;
			}
		}
		if (force_continue)
		{
			continue;
		}
		iter = cur_iters[cur_distance];
		if (*iter == query_node2)
		{
			// add result path and continue;
			
			vector<NODE_TYPE> temp_result_path(c_path_v);
			//vector<NODE_TYPE> temp_result_path(c_path., c_path.end());
			temp_result_path.push_back(*iter);
			result_path.push_back(temp_result_path);
			cur_iters[cur_distance]++;
			continue;
		}
		else if (hot_points.find(*iter) != hot_points.end())
		{
			cur_iters[cur_distance]++;
			continue;
		}

		cur_distance++;
		cur_iters[cur_distance] = adjacency_list[*cur_iters[cur_distance - 1]].begin();

		cur_node = *cur_iters[cur_distance - 1];
		c_path_set.insert(cur_node);
		c_path.push(cur_node);
		c_path_v.push_back(cur_node);
		//cur_iters[cur_distance]++;

	}
	return result_path;
}



void single_direction_baseline_stop_at_hotpoints_bfs_spd(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  k, set<NODE_TYPE> hot_points, spd_index& index)
{
	unordered_set<NODE_TYPE> visited_node;
	DISTANCE_TYPE cur_distance = 0;
	unordered_set<NODE_TYPE> cur_proprogate;
	visited_node.insert(query_node1);
	//cur_proprogate.insert(query_node1);
	for (auto iter = adjacency_list[query_node1].begin(); iter != adjacency_list[query_node1].end(); iter++)
	{
		cur_proprogate.insert(*iter);
	}
	while (cur_distance <= k && !cur_proprogate.empty())
	{
		cur_distance++;
		unordered_set<NODE_TYPE> temp_proprogate;
		NODE_TYPE cur_node;
		for (auto iter = cur_proprogate.begin(); iter != cur_proprogate.end(); iter++)
		{
			cur_node = *iter;
			if (visited_node.find(cur_node) != visited_node.end())//already visited
			{
				continue;
			}
			if (hot_points.find(cur_node) != hot_points.end())// in hot points
			{
				index.hot_push_back(query_node1, cur_node, cur_distance);
			}
			else {
				index.non_hot_push_back(query_node1, cur_node, cur_distance);
			}
			visited_node.insert(cur_node);
			for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
			{
				temp_proprogate.insert(*iter2);
			}
		}
		cur_proprogate = temp_proprogate;
	}
	
}


spd_distance_map single_direction_baseline_stop_at_hotpoints_dfs_spd(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points)
{//original construct index method
	spd_distance_map result;
	stack<NODE_TYPE> c_path;
	vector<NODE_TYPE> c_path_v;
	set<NODE_TYPE> c_path_set;
	int cur_distance = 1;
	vector< vector<NODE_TYPE>::iterator > cur_iters;
	NODE_TYPE src = query_node1;
	c_path.push(src);
	c_path_v.push_back(src);
	vector<int> min_order;
	for (int i = 0; i <= k + 1; i++)
	{
		cur_iters.push_back(adjacency_list[src].begin());
	}
	//set<NODE_TYPE> new_c_path_set(c_path_set);
	c_path_set.insert(src);
	bool dfs_finish = false;
	NODE_TYPE cur_node;
	while (!dfs_finish)
	{// we do not use adjacency_list_double
		NODE_TYPE cur_node, next_node;
		if (cur_distance <= 0)
		{
			break;
		}
		auto iter = cur_iters[cur_distance];
		if (cur_distance == 1)
		{
			if (cur_iters[cur_distance] == adjacency_list[src].end())
			{
				cur_distance--;
				cur_iters[cur_distance]++;
				cur_node = c_path.top();
				c_path.pop();
				c_path_v.pop_back();
				auto iter2 = c_path_set.find(cur_node);
				if (iter2 != c_path_set.end())
				{
					c_path_set.erase(iter2);
				}
				continue;
			}
		}
		else if (cur_iters[cur_distance] == adjacency_list[*cur_iters[cur_distance - 1]].end() || cur_distance > k)
		{
			// go back
			cur_distance--;
			cur_iters[cur_distance]++;
			cur_node = c_path.top();
			c_path.pop();
			c_path_v.pop_back();
			auto iter2 = c_path_set.find(cur_node);
			if (iter2 != c_path_set.end())
			{
				c_path_set.erase(iter2);
			}
			continue;
		}
		cur_node = *iter;
		//normal logic
		bool force_continue = false;
		while (c_path_set.find(*cur_iters[cur_distance]) != c_path_set.end())//*cur_iters[cur_distance] < src || 
		{
			cur_iters[cur_distance]++;
			if (cur_distance == 1)
			{
				if (cur_iters[cur_distance] == adjacency_list[src].end())
				{
					force_continue = true;
					break;
				}
			}
			else if (cur_iters[cur_distance] == adjacency_list[*cur_iters[cur_distance - 1]].end())
			{
				force_continue = true;
				break;
			}
		}
		if (force_continue)
		{
			continue;
		}
		iter = cur_iters[cur_distance];
		if (*iter == query_node2)
		{
			// add result path and continue;

			vector<NODE_TYPE> temp_result_path(c_path_v);
			//vector<NODE_TYPE> temp_result_path(c_path., c_path.end());
			temp_result_path.push_back(*iter);
			result.push_back(*iter, temp_result_path.size()-1);
			cur_iters[cur_distance]++;
			continue;
		}
		else if (hot_points.find(*iter) != hot_points.end())
		{
			cur_iters[cur_distance]++;
			continue;
		}

		cur_distance++;
		cur_iters[cur_distance] = adjacency_list[*cur_iters[cur_distance - 1]].begin();

		cur_node = *cur_iters[cur_distance - 1];
		c_path_set.insert(cur_node);
		c_path.push(cur_node);
		c_path_v.push_back(cur_node);
		//cur_iters[cur_distance]++;

	}
	return result;
}



paths single_direction_baseline_stop_at_hotpoints(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k,set<NODE_TYPE> hot_points)
{//original construct index method
	set<NODE_TYPE> visited_node;
	visited_node.insert(query_node2);
	visited_node.insert(query_node1);
	paths result_path;
	//ofstream outPath;
	//outPath.open("out.txt");
	NODE_TYPE path_count = 0;
	map<NODE_TYPE, vector<NODE_TYPE > >* parents = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];
	set<NODE_TYPE > cur_propagate;
	NODE_TYPE  cur_distance = 1;
	cur_propagate.insert(query_node1);
	while (cur_distance <= k)
	{
		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes

				if (parents[*iter2].find(cur_distance) == parents[*iter2].end()) {
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
				}
				if (*iter2 != query_node2 && (hot_points.find(*iter2) != hot_points.end()))
				{
				}
				else
				{
					if (visited_node.find(*iter2) != visited_node.end())
					{
						continue;
					}
					temp_propagate.insert(*iter2);//insert all the new propagate nodes
					visited_node.insert(*iter2);
				}
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}
	//output all the path
	result_path = find_all_paths_use_parents(adjacency_list,node_num,query_node1,query_node2,k,parents);
	delete[] parents;
	return result_path;
}



//change to use map to construct parents 
paths single_direction_baseline(vector<NODE_TYPE >* adjacency_list, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k)
{
	set<NODE_TYPE> visited_node;
	visited_node.insert(query_node2);
	visited_node.insert(query_node1);
	paths result_path;
	//ofstream outPath;
	//outPath.open("out.txt");
	NODE_TYPE path_count = 0;
	map<NODE_TYPE, vector<NODE_TYPE > >* parents = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];
	set<NODE_TYPE > cur_propagate;
	NODE_TYPE  cur_distance = 1;
	cur_propagate.insert(query_node1);
	while (cur_distance <= k)
	{
		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes
				if (parents[*iter2].find(cur_distance) == parents[*iter2].end()) {
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
				}
				if (visited_node.find(*iter2) != visited_node.end())
				{
					continue;
				}
				temp_propagate.insert(*iter2);//insert all the new propagate nodes
				visited_node.insert(*iter2);
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}
	//output all the path
	result_path = find_all_paths_use_parents(adjacency_list,node_num,query_node1,query_node2,k,parents);
	return result_path;
}


void dfs_ex_node2(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, map<NODE_TYPE, paths> &result, NODE_TYPE cur_node, vector<NODE_TYPE> c_path, set<NODE_TYPE> stop_nodes, int distance, int cur_distance, int& min_stop_disatance, set<NODE_TYPE> c_path_set, NODE_TYPE query_node2)// distance means length of path
{
	//stop condition

	if (stop_nodes.find(cur_node) != stop_nodes.end())
	{
		if (cur_distance < min_stop_disatance && cur_node != query_node2)
		{
			min_stop_disatance = cur_distance;
		}
		if (cur_distance == 0)// for the case the start node is stop node, then return an empty result set
		{
			paths temp_result;
			vector<NODE_TYPE> path;
			path.push_back(cur_node);
			temp_result.push_back(path);
			result.insert(make_pair(cur_node, temp_result));
			return;
		}
		//if meet a stop nodes
		vector<NODE_TYPE> temp_result_path(c_path);
		temp_result_path.push_back(cur_node);
		if (result.find(cur_node) == result.end())
		{//first time to insert result
			paths new_paths;
			new_paths.push_back(temp_result_path);
			result.insert(make_pair(cur_node, new_paths));
		}
		else {
			result.find(cur_node)->second.push_back(temp_result_path);
		}
		return;
	}
	if (cur_distance >= distance)
	{
		return;
	}
	c_path.push_back(cur_node);
	set<NODE_TYPE> new_c_path_set(c_path_set);
	new_c_path_set.insert(cur_node);
	for (auto iter = adjacency_list_double[cur_node].begin(); iter != adjacency_list_double[cur_node].end(); iter++)
	{
		if (c_path_set.find(*iter) != c_path_set.end())// there is a repeat node
		{
			continue;
		}
		dfs_ex_node2(adjacency_list_double, adjacency_list, result, *iter, c_path, stop_nodes, distance, cur_distance + 1, min_stop_disatance, new_c_path_set,query_node2);
		//c_path_set.erase(*iter);
	}

	for (auto iter = adjacency_list[cur_node].begin(); iter != adjacency_list[cur_node].end(); iter++)
	{
		if (c_path_set.find(*iter) != c_path_set.end())// there is a repeat node
		{
			continue;
		}
		dfs_ex_node2(adjacency_list_double, adjacency_list, result, *iter, c_path, stop_nodes, distance, cur_distance + 1, min_stop_disatance, new_c_path_set,query_node2);
		//c_path_set.erase(*iter);
	}
}



void dfs(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, map<NODE_TYPE, paths> &result,NODE_TYPE cur_node, vector<NODE_TYPE> c_path,  set<NODE_TYPE> stop_nodes, int distance, int cur_distance,int& min_stop_disatance, set<NODE_TYPE> c_path_set)// distance means length of path
{
	//stop condition

	if (stop_nodes.find(cur_node) != stop_nodes.end())
	{
		if (cur_distance < min_stop_disatance)
		{
			min_stop_disatance = cur_distance;
		}
		if (cur_distance == 0)// for the case the start node is stop node, then return an empty result set
		{
			paths temp_result;
			vector<NODE_TYPE> path;
			path.push_back(cur_node);
			temp_result.push_back(path);
			result.insert(make_pair(cur_node, temp_result));
			return;
		}
		//if meet a stop nodes
		vector<NODE_TYPE> temp_result_path(c_path);
		temp_result_path.push_back(cur_node);
		if (result.find(cur_node) == result.end())
		{//first time to insert result
			paths new_paths;
			new_paths.push_back(temp_result_path);
			result.insert(make_pair(cur_node, new_paths));
		}
		else {
			result.find(cur_node)->second.push_back(temp_result_path);
		}
		return;
	}
	if (cur_distance >= distance)
	{
		return;
	}
	c_path.push_back(cur_node);
	set<NODE_TYPE> new_c_path_set(c_path_set);
	new_c_path_set.insert(cur_node);
	for (auto iter = adjacency_list_double[cur_node].begin(); iter != adjacency_list_double[cur_node].end(); iter++)
	{
		if (c_path_set.find(*iter) != c_path_set.end())// there is a repeat node
		{
			continue;
		}
		dfs(adjacency_list_double, adjacency_list, result, *iter, c_path, stop_nodes, distance, cur_distance + 1, min_stop_disatance, new_c_path_set);
		//c_path_set.erase(*iter);
	}

	for (auto iter = adjacency_list[cur_node].begin(); iter != adjacency_list[cur_node].end(); iter++)
	{
		if (c_path_set.find(*iter) != c_path_set.end())// there is a repeat node
		{
			continue;
		}
		dfs(adjacency_list_double, adjacency_list, result, *iter, c_path, stop_nodes, distance, cur_distance + 1, min_stop_disatance, new_c_path_set);
		//c_path_set.erase(*iter);
	}
}

path_index construct_hot_point_index_dfs(vector<NODE_TYPE >* adjacency_list_double , vector<NODE_TYPE >* adjacency_list, NODE_TYPE  k, set<NODE_TYPE> hot_points)
{
	// cout << " k is " << k << " in contruct index dfs" << endl;
	path_index index(hot_points);
	// use bfs to find some paths between node A and B
	for (auto start = hot_points.begin(); start != hot_points.end(); start++)
	{
		//// cout << "start to find path " << *start << ":" << *end << endl;
		paths temp;
		//temp = single_direction_baseline(adjacency_list,node_num,*start,*end,k);
		map<NODE_TYPE, paths> result;
		vector<NODE_TYPE> c_path;
		set<NODE_TYPE> c_path_set;
		int min_stop_distance = int(k);
		set<NODE_TYPE> temp_stop_points(hot_points);
		for (auto iter = temp_stop_points.begin(); iter != temp_stop_points.end(); )
		{
			if (*iter == *start) {
				//// cout << *iter << " remove in this " << endl;
				iter = temp_stop_points.erase(iter);
			}
			else {
				//// cout << *iter << " in hot points" << endl;
				iter++;
			}
		}
		//temp_stop_points.insert()
		//pruned_subgraph_unordered_map_double_direction temp_pruned;
		//temp_pruned.meet_nodes = temp_stop_points;
		//result = temp_pruned.dfs_without_recursion_all_meetpoints()
		dfs(adjacency_list_double,adjacency_list,result,*start,c_path, temp_stop_points,k,0, min_stop_distance, c_path_set);
		// cout << *start << " :  in construct bfs "<<endl;
		for (auto iter = result.begin(); iter != result.end(); iter++)
		{
			index.push_back(iter->second);
		}
		
	}
	return index;
	// finish index construct

}

path_index construct_hot_point_index_dfs_using_new_algorithm(unordered_map<NODE_TYPE, set<NODE_TYPE>> &induced_subgraph, NODE_TYPE  k, set<NODE_TYPE> hot_points, NODE_TYPE node_num)
{
	// cout << " k is " << k << " in contruct index dfs" << endl;
	path_index index(hot_points);
	// use bfs to find some paths between node A and B
	for (auto start = hot_points.begin(); start != hot_points.end(); start++)
	{
		//// cout << "start to find path " << *start << ":" << *end << endl;
		paths temp;
		//temp = single_direction_baseline(adjacency_list,node_num,*start,*end,k);
		map<NODE_TYPE, paths> result;
		vector<NODE_TYPE> c_path;
		set<NODE_TYPE> c_path_set;
		int min_stop_distance = int(k);
		set<NODE_TYPE> temp_stop_points(hot_points);
		for (auto iter = temp_stop_points.begin(); iter != temp_stop_points.end(); )
		{
			if (*iter == *start) {
				//// cout << *iter << " remove in this " << endl;
				iter = temp_stop_points.erase(iter);
			}
			else {
				//// cout << *iter << " in hot points" << endl;
				iter++;
			}
		}
		pruned_subgraph_unordered_map_double_direction temp_pruned;
		temp_pruned.meet_nodes = temp_stop_points;
		result = temp_pruned.dfs_without_recursion_all_meetpoints(induced_subgraph, *start, NULL_NODE, k, node_num);
		//dfs(adjacency_list_double, adjacency_list, result, *start, c_path, temp_stop_points, k, 0, min_stop_distance, c_path_set);
		// cout << *start << " :  in construct bfs " << endl;
		for (auto iter = result.begin(); iter != result.end(); iter++)
		{
			index.push_back(iter->second);
		}

	}
	return index;
	// finish index construct

}



path_index construct_hot_point_index(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  k, set<NODE_TYPE> hot_points)
{
	// cout << " k is " << k << " in contruct index " << endl;
	path_index index(hot_points);
	// use bfs to find some paths between node A and B
	for (auto start = hot_points.begin(); start != hot_points.end(); start++)
	{
		for (auto end = hot_points.begin(); end != hot_points.end(); end++)
		{
			if (*start == *end)
			{
				continue;
			}
			//// cout << "start to find path " << *start << ":" << *end << endl;
			paths temp;;
			//temp = single_direction_baseline(adjacency_list,node_num,*start,*end,k);
			temp = single_direction_baseline_stop_at_hotpoints_dfs(adjacency_list, node_num, *start, *end, k,hot_points);
			//// cout << *start << " : " << *end << " in bfs "<<endl;
	
			index.push_back(temp);
		}
	}
	return index;
	// finish index construct

}


set<NODE_TYPE> find_hot_points_top_t(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  k, int t)
{
	vector<hot_degree> hot_points;
	//double threhold = 0.05;//this means we find top 5% nodes in degree rank as hot-points(this part can change to k-core value, truess value or centrality)
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		NODE_TYPE temp_degree = adjacency_list[i].size() + adjacency_list_double[i].size() +adjacency_list_reverse[i].size();// sum in and out degree
		if (temp_degree != 0)
		{
			hot_degree hd1(i, temp_degree);
			hot_points.push_back(hd1);
		}
	}
	stable_sort(hot_points.begin(), hot_points.end(), less<hot_degree>());
	//NODE_TYPE end = (NODE_TYPE)(hot_points.size() * threshold);
	//// cout << " total nodes is " << hot_points.size() << endl;
	//// cout << " we choose " << end << endl;
	set<NODE_TYPE> result;
	for (NODE_TYPE i = 0; i < t; i++)
	{
		result.insert(hot_points[i].node);
		//// cout << hot_points[i].node << endl;
	}
	// cout << "hot points : " << result.size();
	//for (auto iter = result.begin(); iter != result.end(); iter++)
	//{
	//	// cout << *iter << " ";
	//}
	// cout << endl;
	return result;
}

set<NODE_TYPE> find_hot_points(vector<NODE_TYPE >* adjacency_list_double,vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  k, double threshold)
{
	vector<hot_degree> hot_points;
	//double threhold = 0.05;//this means we find top 5% nodes in degree rank as hot-points(this part can change to k-core value, truess value or centrality)
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		NODE_TYPE temp_degree = adjacency_list[i].size() + adjacency_list_double[i].size();// +adjacency_list_reverse[i].size();// sum in and out degree
		if (temp_degree != 0)
		{
			hot_degree hd1(i, temp_degree);
			hot_points.push_back(hd1);
		}
	}
	stable_sort(hot_points.begin(), hot_points.end(), less<hot_degree>());
	NODE_TYPE end = (NODE_TYPE)(hot_points.size() * threshold);
	//// cout << " total nodes is " << hot_points.size() << endl;
	//// cout << " we choose " << end << endl;
	set<NODE_TYPE> result;
	for (NODE_TYPE i = 0; i < end; i++)
	{
		result.insert(hot_points[i].node);
		//// cout << hot_points[i].node << endl;
	}
	// cout << "hot points : " << result.size();
	//for (auto iter = result.begin(); iter != result.end(); iter++)
	//{
	//	// cout << *iter << " ";
	//}
	// cout << endl;
	return result;
}




void printpath(vector<NODE_TYPE>& path)
{
	int size = path.size();
	for (int i = 0; i < size; i++)
	{
		cout << path[i] << " ";
	}
	cout << endl;
}

// utility function to check if current
// vertex is already present in path
int isNotVisited(NODE_TYPE x, vector<NODE_TYPE>& path)
{
	int size = path.size();
	for (int i = 0; i < size; i++)
		if (path[i] == x)
			return 0;
	return 1;
}

// utility function for finding paths in graph
// from source to destination

//rewrite bfs for find all paths
two_nodes_path_index findpaths(vector<NODE_TYPE>* g, NODE_TYPE src, NODE_TYPE dst, int k)
{
	// create a queue which stores
	// the paths
	queue<vector<NODE_TYPE> > q;
	two_nodes_path_index result(src, dst);
	//vector<vector<NODE_TYPE> > result;
	// path vector to store the current path
	vector<NODE_TYPE> path;
	path.push_back(src);
	q.push(path);
	while (!q.empty()) {
		path = q.front();
		q.pop();
		
		int last = path[path.size() - 1];

		// if last vertex is the desired destination
		// then print the path
		if (last == dst) {
			//printpath(path);
			result.push_back(path);
			// path remove duplicate
		}
		if (path.size() >= k+1)// a path with length k will have k+1 nodes
		{
			continue;
		}
		// traverse to all the nodes connected to 
		// current vertex and push new path to queue
		for (int i = 0; i < g[last].size(); i++) {
			if (isNotVisited(g[last][i], path)) {
				vector<NODE_TYPE> newpath(path);
				newpath.push_back(g[last][i]);
				q.push(newpath);
			}
		}
	}
	//result.output();
	return result;
}

void test_baseline_and_index_performance()
{
	set<pair<NODE_TYPE, NODE_TYPE > > a;
	NODE_TYPE query1, query2, k;
	//query1 = 1;
	//query2 = 498948;
	//query1 = 222074;
	//query2 = 458358;
	//test1

	//query1 = 522644;
	//query2 = 498948;
	//k = 6;

	//query1 = 233;
	//query2 = 456;
	//k = 6;

	//query1 = 222074;
	//query2 = 458358;
	//k = 4;

	query1 = 1;
	query2 = 6;
	k = 5;

	NODE_TYPE  node_num = 1007518272;//548552;//1234944765;
								  //NODE_TYPE node_num = 91306;
	vector<NODE_TYPE >* adjacency_list = new vector<NODE_TYPE >[node_num + 1];
	vector<NODE_TYPE >* adjacency_list_reverse = new vector<NODE_TYPE >[node_num + 1];
	vector<NODE_TYPE >* adjacency_list_double = new vector<NODE_TYPE>[node_num + 1];//record the double directed adjacency for the undirected graph load 
																					//load_graph_data("Amazon.txt", adjacency_list, adjacency_list_reverse);
	//load_undirected_graph_data("test_case.txt", adjacency_list_double);
	load_ali_undirected_data(adjacency_list_double, adjacency_list, adjacency_list_reverse);
	clock_t index_start, index_end;
	index_start = clock();
	set<NODE_TYPE> hot_points = find_hot_points(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, k,0.05);
	// cout << hot_points.size() << " is the number of hot points " << endl;
	//hot_points.insert(1);
	//hot_points.insert(6);
	path_index index;
	//index = construct_hot_point_index(adjacency_list, adjacency_list_reverse, node_num, k, hot_points);
	index = construct_hot_point_index_dfs(adjacency_list_double, adjacency_list, k, hot_points);
	index.output();
	// cout << " end of index " << endl;
	index_end = clock();


	//
	int table = 1;
	NODE_TYPE x, y;
	string table_name;
	ifstream inEdges;
	ofstream performance;
	performance.open("performance.txt");
	char buffer[BUFFER_LENTH];
	NODE_TYPE query_index = 0;
	for (table = 0; table <= 19; table++)
	{
		table_name = "./data/dwd_rsm_dmt_graph_data20171228_04_new_dd_idx/dwd_rsm_dmt_graph_data20171228_04_new_dd_idx_" + to_string(table);
		// cout << table_name << " is our talbe name" << endl;
		//for each table, input all the key-value pairs. Store and output all ordered key-value pairs.
		inEdges.open(table_name);
		if (!inEdges.is_open())
		{
			// cout << "Error opening file : " << table_name << endl;
			continue;
		}
		while (!inEdges.eof())
		{
			query_index++;
			inEdges.getline(buffer, BUFFER_LENTH);
			extractEdges(buffer, x, y);
			clock_t start, start1, end, end2;
			paths single;
			map<NODE_TYPE, paths> dfs_result;
			vector<NODE_TYPE> c_path;
			set<NODE_TYPE> c_path_set;
			set<NODE_TYPE> stop_points;
			int min_stop_distance = k;
			query1 = y;
			query2 = x;
			stop_points.insert(query2);
			//single = single_direction_baseline(adjacency_list, node_num, query1, query2, k);
			start = clock();
			dfs(adjacency_list_double, adjacency_list, dfs_result, query1, c_path, stop_points, k, 0, min_stop_distance, c_path_set);
			end = clock();
			if (dfs_result.find(query2) == dfs_result.end()) {
				paths empty_path;
				empty_path.write_to_file("./result/single" + to_string(query_index) + ".txt");
			}
			else {
				dfs_result.find(query2)->second.write_to_file("./result/single" + to_string(query_index) + ".txt");
			}
			paths h_result2;
			start1 = clock();
			h_result2 = find_all_paths_with_hotpoints_dfs(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, query1, query2, k, hot_points, index);
			end2 = clock();
			h_result2.write_to_file("./result/index_result_step"+to_string(query_index)+".txt");
			//adjacency_list_double[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list
												  //adjacency_list_reverse[y].push_back(x);
			performance << x <<"\t:\t" <<y<<"\t"<<double(end - start) / CLOCKS_PER_SEC << "\t:\t" << double(end2 - start1) / CLOCKS_PER_SEC << endl;
		}
		inEdges.close();
	}
	performance.close();

}



paths double_direction_baseline_paths(vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k)
{
	paths results;
	map<NODE_TYPE, paths > first_path_map;
	map<NODE_TYPE, paths > second_path_map;
	NODE_TYPE path_count = 0;
	bool* is_board_nodes = new bool[node_num + 1];
	for (int i = 0; i <= node_num; i++)
	{
		is_board_nodes[i] = false;
	}
	map<NODE_TYPE, NODE_TYPE> path_count_map;
	map<NODE_TYPE, NODE_TYPE> path_count_map_second;
	map<NODE_TYPE, vector<NODE_TYPE > >* parents = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];//record distance and last nodes
	set<NODE_TYPE > cur_propagate;
	NODE_TYPE  cur_distance = 1;
	cur_propagate.insert(query_node1);
	while (cur_distance <= int(k / 2))
	{
		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			//// cout << "cur_node is " << cur_node << endl;
			for (auto iter2 = adjacency_list[cur_node].begin(); iter2 != adjacency_list[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes
			 //// cout << "iter2 " << *iter2 << endl;

				if (parents[*iter2].find(cur_distance) == parents[*iter2].end()) {
					//if there is no vector for this distance, pair is distance, node 
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
					//parents_meetup[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
													  //iter3 = parents_meetup[*iter2].find(cur_distance);
													  //iter3->second.push_back(cur_node);//insert node distance value
				}
				//		auto iter4 = parents[*iter2].find(cur_distance);
				if (cur_distance == int(k / 2))//record all the board nodes
				{
					is_board_nodes[*iter2] = true;
				}
				temp_propagate.insert(*iter2);//insert all the new propagate nodes
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}

	cur_propagate.clear();
	cur_propagate.insert(query_node2);//search reversely
	cur_distance = 1;
	map<NODE_TYPE, vector<NODE_TYPE > >* parents_query2 = new map<NODE_TYPE, vector<NODE_TYPE > >[node_num + 1];//record distance and last nodes
	stack<pair<NODE_TYPE, NODE_TYPE > > short_path;
	if (!parents[query_node2].empty())
	{
		//output all the short path and remove all the nodes on these path
		for (auto iter = parents[query_node2].begin(); iter != parents[query_node2].end(); iter++)
		{
			for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
			{
				short_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter->first));
				//parents_meetup[query_node2].erase(iter->first);

			}
		}
	}
	vector<NODE_TYPE> output_path;
	output_path.push_back(query_node2);
	if (first_path_map.find(query_node2) == first_path_map.end())
	{
		paths a;
		first_path_map.insert(pair<NODE_TYPE, paths >(query_node2, a));
		//a.push_back("");
		second_path_map.insert(pair<NODE_TYPE, paths >(query_node2, a));
	}


	//// cout << "short path is " << endl;
	while (!short_path.empty())
	{
		NODE_TYPE  temp_node = short_path.top().first;
		NODE_TYPE  distance = short_path.top().second;
		if (temp_node == query_node1)
		{
			//output_path
			vector<NODE_TYPE> temp_path_vector;
			for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
			{
				//// cout << *iter4 << " ";
				//out_path << *iter4 << " ";
				temp_path_vector.push_back(*iter4);
			}
			//// cout << query_node1 << endl;
			//out_path << query_node1 << endl;
			temp_path_vector.push_back(query_node1);
			first_path_map.find(query_node2)->second.push_back(temp_path_vector);
			vector<NODE_TYPE> temp_result(temp_path_vector);
			std::reverse(temp_result.begin(), temp_result.end());
			results.push_back(temp_result);
			path_count++;
			//output_path.pop_back();
			short_path.pop();
			continue;
		}
		else if (temp_node == output_path.back())
		{
			output_path.pop_back();
			short_path.pop();
			continue;
		}
		else {
			output_path.push_back(temp_node);
		}
		//path.pop();
		/*if (distance == k)
		{
		// cout << endl;
		}*/
		//// cout << temp_node<<":" << distance << " ";
		distance--;
		auto iter = parents[temp_node].find(distance);
		/*if (distance == 0 || iter == parents[temp_node].end())
		{
		// cout << " end of path" << endl;
		continue;
		}*/
		auto iter2 = iter->second;
		for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
		{
			short_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
		}
	}
	set<NODE_TYPE > another_path;
	is_board_nodes[query_node2] = false;// because query_node2 is outputed in the short path section
	is_board_nodes[query_node1] = false;
	while (cur_distance <= k - int(k / 2))// meet up reversely 
	{
		//start to output all the meetup paths

		set<NODE_TYPE > temp_propagate;
		for (auto iter = cur_propagate.begin(); iter != cur_propagate.end(); iter++)//traverse all the cur_propagate set nodes
		{//insert <node,distance> pair first
			NODE_TYPE  cur_node = *iter;
			//// cout << "cur_node is " << cur_node << endl;
			for (auto iter2 = adjacency_list_reverse[cur_node].begin(); iter2 != adjacency_list_reverse[cur_node].end(); iter2++)
			{//traverse all the adjacent nodes

				if (parents_query2[*iter2].find(cur_distance) == parents_query2[*iter2].end()) {
					vector<NODE_TYPE > a;
					a.push_back(cur_node);
					parents_query2[*iter2].insert(pair<NODE_TYPE, vector<NODE_TYPE > >(cur_distance, a));
				}
				else {
					auto iter3 = parents_query2[*iter2].find(cur_distance);
					iter3->second.push_back(cur_node);//insert node distance value
				}
				//		auto iter4 = parents[*iter2].find(cur_distance);

				temp_propagate.insert(*iter2);//insert all the new propagate nodes
											  //judge whether there are meet-up paths at cur_node
				if (is_board_nodes[*iter2] == true)
				{
					if (first_path_map.find(*iter2) == first_path_map.end())
					{
						paths a;
						first_path_map.insert(pair<NODE_TYPE, paths >(*iter2, a));
						second_path_map.insert(pair<NODE_TYPE, paths >(*iter2, a));
					}
					//output the first half path
					is_board_nodes[*iter2] = false;
					//// cout << "meetup node is " << *iter2 << endl;
					//out_path << "meetup node is " << *iter2 << endl;
					path_count_map.insert(pair<NODE_TYPE, NODE_TYPE>(*iter2, 0));
					path_count_map_second.insert(pair<NODE_TYPE, NODE_TYPE>(*iter2, 0));
					//// cout << endl << "first half path: " << endl;
					//out_path << endl << "first half path: " << endl;
					stack<pair<NODE_TYPE, NODE_TYPE > > temp_path;
					while (!temp_path.empty())
					{
						temp_path.pop();
					}
					//output all the short path and remove all the nodes on these path
					auto iter_temp_board = parents[*iter2].find(int(k / 2));
					if (iter_temp_board != parents[*iter2].end())
					{
						for (auto iter_board = iter_temp_board->second.begin(); iter_board != iter_temp_board->second.end(); iter_board++)
						{
							temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter_board, int(k / 2)));
						}
					}

					//temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2,int(k/2)+1 ) ); // distacne node pair
					/*for (auto iter_meet = parents_meetup[cur_node].begin(); iter_meet != parents_meetup[cur_node].end(); iter++)
					{
					for (auto iter2 = iter_meet->second.begin(); iter2 != iter_meet->second.end(); iter2++)
					{
					temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter_meet->first));
					}
					}*/
					//parents_meetup[*iter2].clear();
					auto temp_iter = path_count_map.find(*iter2);
					vector<NODE_TYPE> output_path;
					//output_path.push_back(*iter2);// to be confirm
					while (!temp_path.empty())
					{
						NODE_TYPE  temp_node = temp_path.top().first;
						NODE_TYPE  distance = temp_path.top().second;
						if (temp_node == query_node1)
						{
							//output_path
							vector<NODE_TYPE> temp_path_vector;
							//if (!circle_detection(output_path))
							{
								for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
								{
									//// cout << *iter4 << " ";
									//out_path << *iter4 << " ";
									temp_path_vector.push_back(*iter4);
								}
								//// cout << query_node1 << endl;
								//out_path << query_node1 << endl;
								temp_path_vector.push_back(query_node1);
								first_path_map.find(*iter2)->second.push_back(temp_path_vector);
								temp_iter->second += 1;
								//output_path.pop_back();
							}
							temp_path.pop();
							continue;
						}
						else if (output_path.size() == 0)
						{
							output_path.push_back(temp_node);
						}
						else if (temp_node == output_path.back())
						{
							output_path.pop_back();
							temp_path.pop();
							continue;
						}
						else {
							output_path.push_back(temp_node);
						}
						//temp_path.pop();
						//// cout << temp_node << " ";
						distance--;
						if (parents[temp_node].find(distance) == parents[temp_node].end())
						{
							while (!temp_path.empty())
							{
								temp_path.pop();
							}
							break;
						}
						auto iter_parent = parents[temp_node].find(distance);
						//auto iter3 = parents_meetup[temp_node].find(distance);
						/*if (distance == 0 || iter == parents[temp_node].end())
						{
						// cout << " end of path" << endl;
						continue;
						}*/
						auto iter2 = iter_parent->second;
						//iter3->second.clear();
						//parents_meetup[temp_node].erase(distance);
						for (auto iter2 = iter_parent->second.begin(); iter2 != iter_parent->second.end(); iter2++)
						{
							temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
						}
					}

					// output another half path
					another_path.insert(*iter2);

				}
			}
		}
		cur_propagate = temp_propagate;
		cur_distance++;
	}

	//// cout << "another half path" << endl;
	//out_path << "another half path" << endl;
	stack<pair<NODE_TYPE, NODE_TYPE > > temp_path;
	for (auto iter_set = another_path.begin(); iter_set != another_path.end(); iter_set++)
	{
		NODE_TYPE  cur_node = *iter_set;
		//// cout << "meet node is " << cur_node << endl;
		//out_path << "meet node is " << cur_node << endl;
		if (!parents_query2[cur_node].empty())
		{
			//output all the short path and remove all the nodes on these path
			for (auto iter = parents_query2[cur_node].begin(); iter != parents_query2[cur_node].end(); iter++)
			{
				for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
				{
					temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, iter->first));
				}
			}
		}
		auto temp_iter2 = path_count_map_second.find(cur_node);
		vector<NODE_TYPE> output_path;
		output_path.push_back(cur_node);
		while (!temp_path.empty())
		{
			NODE_TYPE  temp_node = temp_path.top().first;
			NODE_TYPE  distance = temp_path.top().second;
			//temp_path.pop();
			//// cout << temp_node << " ";
			if (temp_node == query_node2)
			{
				vector<NODE_TYPE> temp_path_vector;
				//output_path
				for (auto iter4 = output_path.begin(); iter4 != output_path.end(); iter4++)
				{
					//// cout << *iter4 << " ";
					//out_path << *iter4 << " ";
					temp_path_vector.push_back(*iter4);
				}
				//// cout << query_node2 << endl;
				//out_path << query_node2 << endl;
				temp_path_vector.push_back(query_node2);
				second_path_map.find(cur_node)->second.push_back(temp_path_vector);
				temp_iter2->second += 1;
				//output_path.pop_back();
				temp_path.pop();
				continue;
			}
			else if (temp_node == output_path.back())
			{
				output_path.pop_back();
				temp_path.pop();
				continue;
			}
			else {
				output_path.push_back(temp_node);
			}
			distance--;
			auto iter = parents_query2[temp_node].find(distance);
			//auto iter3 = parents_meetup[temp_node].find(distance);
			/*if (distance == 0 || iter == parents_query2[temp_node].end())
			{
			// cout << " end of path" << endl;
			continue;
			}*/
			auto iter2 = iter->second;
			//iter3->second.clear();
			for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++)
			{
				temp_path.push(pair<NODE_TYPE, NODE_TYPE >(*iter2, distance));
			}
		}
	}
	for (auto iter_count = path_count_map.begin(); iter_count != path_count_map.end(); iter_count++)
	{
		auto iter_count2 = path_count_map_second.find(iter_count->first);
		//// cout << iter_count->first << " is meet point" << iter_count2->second << ":" << iter_count->second << " in total" << endl;
		path_count += iter_count2->second * iter_count->second;
	}

	//// cout << "all the path" << endl;
	//out_path << "all the path" << endl;
	for (auto iter = first_path_map.begin(); iter != first_path_map.end(); iter++)
	{
		auto iter2 = second_path_map.find(iter->first);
		for (auto iter3 = iter->second.path.begin(); iter3 != iter->second.path.end(); iter3++)
		{
			for (auto iter4 = iter2->second.path.begin(); iter4 != iter2->second.path.end(); iter4++)
			{
				//if (!circle_detection_string(reverse_path(*iter3) + *iter4))
				//{
					//// cout << reverse_path(*iter3) << *iter4 << endl;
					//vector<NODE_TYPE> temp_result;
					//vector<NODE_TYPE> temp2 = reverse_path(*iter3);
					//temp_result.insert(temp_result.end()), 
					//result.push_back(reverse_path(*iter3) + (*iter4));
				vector<NODE_TYPE> temp(*iter3);
				std::reverse(temp.begin(), temp.end());
				temp.insert(temp.end(), iter4->begin(), iter4->end());
				results.push_back(temp);
					//out_path << reverse_path(*iter3) << *iter4 << " " << endl;
				//}
			}
		}
	}


	//out_path.close();
	//// cout << path_count << " path in double direction method" << endl;
	
	results.drop_path_length_more_than_k(k);
	results.drop_path_with_repeat_node();
	results.sort_by_string_order();
	results.drop_repeat_path();
	delete[] is_board_nodes;
	delete[] parents;
	delete[] parents_query2;
	return results;

}
