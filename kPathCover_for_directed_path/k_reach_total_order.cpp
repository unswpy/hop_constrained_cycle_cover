#include "k_reach_total_order.h"
#include "BFS_prune.h"
void set_new_total_order(NODE_TYPE node_num, vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, NODE_TYPE* map_table)
{
	vector<NODE_TYPE>* adjacency_list2 = new vector<NODE_TYPE>[node_num + 1];
	vector<NODE_TYPE>* adjacency_list_reverse2 = new vector<NODE_TYPE>[node_num + 1];
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		adjacency_list2[map_table[i]] = adjacency_list[i];
		for (auto iter = adjacency_list[i].begin(); iter != adjacency_list[i].end(); iter++)
		{
			*iter = map_table[*iter];
		}
	}
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		adjacency_list_reverse2[map_table[i]] = adjacency_list_reverse[i];
		for (auto iter = adjacency_list_reverse[i].begin(); iter != adjacency_list_reverse[i].end(); iter++)
		{
			*iter = map_table[*iter];
		}
	}
	delete[] adjacency_list;
	delete[] adjacency_list_reverse;
	adjacency_list = adjacency_list2;
	adjacency_list_reverse = adjacency_list_reverse2;
}

unordered_map<NODE_TYPE, DISTANCE_TYPE> BFS_from_t(NODE_TYPE t, vector<NODE_TYPE>* adjacency_list_reverse, int k, set<NODE_TYPE> block_nodes)
{//return k-reach nodes from source node s
	unordered_map<NODE_TYPE, DISTANCE_TYPE> result;
	set<NODE_TYPE> visited;
	set<NODE_TYPE> cur_propagate;
	int cur_dis = 0;
	result.insert(std::make_pair(t, 0));
	visited.insert(t);
	cur_propagate.insert(t);
	for (cur_dis = 0; cur_dis < k; cur_dis++)
	{
		set<NODE_TYPE> temp_proprogate;
		for (auto iter1 = cur_propagate.begin(); iter1 != cur_propagate.end(); iter1++)
		{
			for (auto iter2 = adjacency_list_reverse[*iter1].begin(); iter2 != adjacency_list_reverse[*iter1].end(); iter2++)
			{
				if (block_nodes.find(*iter2) == block_nodes.end() && visited.find(*iter2) == visited.end() )// if *iter2 is not in block set
				{
					result.insert(std::make_pair(*iter2, cur_dis));
					temp_proprogate.insert(*iter2);
					visited.insert(*iter2);
				}
			}
		}
		cur_propagate = temp_proprogate;
	}
	return result;
}

void list_paths_theoryG(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, NODE_TYPE cur_node, NODE_TYPE query_node1, NODE_TYPE query_node2, int k, paths &result, set<NODE_TYPE> c_path_set, int cur_distance, vector<NODE_TYPE> c_path)
{
	if (cur_node == query_node2)
	{

		//if meet a stop nodes
		vector<NODE_TYPE> temp_result_path(c_path);
		temp_result_path.push_back(cur_node);
		result.push_back(temp_result_path);
		return;
	}
	if (cur_distance >= k)
	{
		return;
	}
	c_path.push_back(cur_node);
	set<NODE_TYPE> new_c_path_set(c_path_set);
	new_c_path_set.insert(cur_node);
	unordered_map<NODE_TYPE, DISTANCE_TYPE> dst_t;
	dst_t = BFS_from_t(query_node2, adjacency_list_reverse, k, new_c_path_set);
	for (auto iter = adjacency_list[cur_node].begin(); iter != adjacency_list[cur_node].end(); iter++)
	{
		if (c_path_set.find(*iter) != c_path_set.end())// there is a repeat node or *iter cannot reach dst
		{
			continue;
		}
		auto iter_dst = dst_t.find(*iter);
		if (iter_dst == dst_t.end())
		{
			continue;
		}
		else if (iter_dst->second + cur_distance > k)
		{
			continue;
		}
		list_paths_theoryG(adjacency_list, adjacency_list_reverse, *iter, query_node1, query_node2, k, result, new_c_path_set, cur_distance + 1, c_path);
		//c_path_set.erase(*iter);
	}
}

void list_paths_theoryG_plus(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, NODE_TYPE cur_node, NODE_TYPE query_node1, NODE_TYPE query_node2, int k, paths &result, set<NODE_TYPE> c_path_set, int cur_distance, vector<NODE_TYPE> c_path)
{//at initial stage, we need to find a shortest path from s to t
	if (cur_node == query_node2)
	{
		//if meet a stop nodes
		vector<NODE_TYPE> temp_result_path(c_path);
		temp_result_path.push_back(cur_node);
		result.push_back(temp_result_path);
		return;
	}
	if (cur_distance >= k)
	{
		return;
	}
	//find shortest path from querynode1 to querynode2

	c_path.push_back(cur_node);
	set<NODE_TYPE> new_c_path_set(c_path_set);
	new_c_path_set.insert(cur_node);
	unordered_map<NODE_TYPE, DISTANCE_TYPE> dst_t;
	dst_t = BFS_from_t(query_node2, adjacency_list_reverse, k, new_c_path_set);
	set<NODE_TYPE> good_neighbors;
	NODE_TYPE temp_cur_node;
	while (true)
	{
		if (cur_distance >= k)
		{
			break;
		}
		good_neighbors.clear();
		for (auto iter = adjacency_list[cur_node].begin(); iter != adjacency_list[cur_node].end(); iter++)
		{
			if (c_path_set.find(*iter) != c_path_set.end())// there is a repeat node or *iter cannot reach dst
			{
				continue;
			}
			auto iter_dst = dst_t.find(*iter);
			if (iter_dst == dst_t.end())
			{
				continue;
			}
			else if (iter_dst->second + cur_distance > k)
			{
				continue;
			}
			if (*iter == query_node2)
			{
				list_paths_theoryG_plus(adjacency_list, adjacency_list_reverse, *iter, query_node1, query_node2, k, result, new_c_path_set, cur_distance + 1, c_path);
			}
			else
			{
				good_neighbors.insert(*iter);
				temp_cur_node = *iter;
			}
			//c_path_set.erase(*iter);
		}
		if (good_neighbors.size() > 1)
		{
			for (auto iter_neighbor = good_neighbors.begin(); iter_neighbor != good_neighbors.end(); iter_neighbor++)
			{
				list_paths_theoryG_plus(adjacency_list, adjacency_list_reverse, *iter_neighbor, query_node1, query_node2, k, result, new_c_path_set, cur_distance + 1, c_path);
			}
			break;//jump out the loop
		}
		else if (good_neighbors.empty())
		{
			break;
		}
		else {
			c_path_set.insert(cur_node);
			cur_node = temp_cur_node;
			cur_distance++;
			c_path.push_back(cur_node);
			new_c_path_set.insert(cur_node);
			
		}
		//c_path.push_back(cur_node);
		
	}

}

void dfs_find_k_paths(vector<NODE_TYPE>* adjacency_list, NODE_TYPE cur_node, NODE_TYPE query_node1, NODE_TYPE query_node2, int k,  paths &result, set<NODE_TYPE> c_path_set, int cur_distance, vector<NODE_TYPE> c_path)
{
	if (cur_node == query_node2)
	{

		//if meet a stop nodes
		vector<NODE_TYPE> temp_result_path(c_path);
		temp_result_path.push_back(cur_node);
		result.push_back(temp_result_path);
		return;
	}
	if (cur_distance >= k)
	{
		return;
	}
	c_path.push_back(cur_node);
	set<NODE_TYPE> new_c_path_set(c_path_set);
	new_c_path_set.insert(cur_node);

	for (auto iter = adjacency_list[cur_node].begin(); iter != adjacency_list[cur_node].end(); iter++)
	{
		if (c_path_set.find(*iter) != c_path_set.end())// there is a repeat node or *iter cannot reach dst
		{
			continue;
		}
		
		dfs_find_k_paths(adjacency_list, *iter,query_node1,query_node2,k,result, new_c_path_set,cur_distance+1,c_path);
		//c_path_set.erase(*iter);
	}
}



NODE_TYPE* get_reverse_map_table(NODE_TYPE node_num, NODE_TYPE* map_table)
{// assume max index of map_table and reverse_map_table is the same
	NODE_TYPE* reverse_map_table = new NODE_TYPE[node_num+1];
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		reverse_map_table[i] = i;
	}
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		reverse_map_table[map_table[i]] = i;
	}
	return reverse_map_table;
}

void construct_k_total_index(NODE_TYPE node_num, vector<map_distance_node_pair> &in_label, vector<map_distance_node_pair> &out_label, vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, vector<NODE_TYPE>* adjacency_list_double, vector<int> in_degrees, vector<int> out_degrees, int k, vector<NODE_TYPE>& map_node_to_total_order, int& cur_low_total_order, set<NODE_TYPE> label_nodes)
{
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		in_label[i].insert(i, 0);
		out_label[i].insert(i, 0);
	}
	for (auto iter = label_nodes.begin(); iter != label_nodes.end(); iter ++)// label nodes may be not in order
	{
		NODE_TYPE i = *iter;
		if (i % 10000 == 0)
		{
			// cout << i << " in " << endl;
		}
		vector<NODE_TYPE> c_path;
		set<NODE_TYPE> c_path_set;
		in_label[i].insert(i, 0);
		if (in_degrees[i] == 0)
		{
			continue;
		}
		construct_k_total_index_dfs_without_recursion(adjacency_list_double, adjacency_list_reverse, i, k, out_label, in_label, node_num,map_node_to_total_order,cur_low_total_order);
		//construct_k_total_index_dfs(adjacency_list_double, adjacency_list_reverse, i, c_path, k, i, 0, c_path_set, out_label, in_label, -1);
	}
	//for (NODE_TYPE i = 0; i < cur_low_total_order; i++)
	for (auto iter = label_nodes.begin(); iter != label_nodes.end(); iter++)
	{
		NODE_TYPE i = *iter;
		if (i == 531483)
		{
			// cout << i << " out " << endl;
		}
		vector<NODE_TYPE> c_path;
		set<NODE_TYPE> c_path_set;
		out_label[i].insert(i, 0);

		if (out_degrees[i] == 0)
		{
			continue;
		}
		construct_k_total_index_dfs_without_recursion(adjacency_list_double, adjacency_list, i, k, in_label, out_label, node_num,map_node_to_total_order,cur_low_total_order);
		//construct_k_total_index_dfs(adjacency_list_double, adjacency_list, i, c_path, k, i, 0, c_path_set, in_label, out_label, -1);
	}



}


// need to debug
void construct_k_total_index_dfs_without_recursion(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, NODE_TYPE src,int k,vector<map_distance_node_pair> &in_label, vector<map_distance_node_pair> &out_label, int node_num, vector<NODE_TYPE>& map_node_to_total_order, int& cur_low_total_order)// distance means length of path
{
	stack<NODE_TYPE> c_path;
	set<NODE_TYPE> c_path_set;
	int cur_distance = 1;
	vector< vector<NODE_TYPE>::iterator > cur_iters;

	c_path.push(src);
	vector<int> min_order;
	for (int i = 0; i <= k+1; i++)
	{
		min_order.push_back(map_node_to_total_order[src]);

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
				min_order[cur_distance] = map_node_to_total_order[src];
				cur_distance--;
				cur_iters[cur_distance]++;
				cur_node = c_path.top();
				c_path.pop();
				auto iter2 = c_path_set.find(cur_node);
				if (iter2 != c_path_set.end())
				{
					c_path_set.erase(iter2);
				}
				continue;
			}
		}
		else if (cur_iters[cur_distance] == adjacency_list[ *cur_iters[cur_distance-1]].end() || cur_distance > k)
		{
			// go back
			min_order[cur_distance-1] = map_node_to_total_order[src];
			cur_distance--;
			cur_iters[cur_distance]++;
			cur_node = c_path.top();
			c_path.pop();
			auto iter2 = c_path_set.find(cur_node);
			if(iter2 != c_path_set.end())
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
		if(cur_distance >= 1 && min_order[cur_distance-1] < min_order[cur_distance])
		{
			min_order[cur_distance] = min_order[cur_distance-1];
		}
		if (map_node_to_total_order[*iter] < min_order[cur_distance])
		{
			min_order[cur_distance] = map_node_to_total_order[*iter];
			if (!out_label[src].intersection(in_label[*iter], k))//if doesnot exist common node
			{
				out_label[src].insert(*iter, cur_distance);
				//in_label[*iter].insert(src, cur_distance);
			}
			//if (!out_label[src].intersection(in_label[*iter], k))//if doesnot exist common node
			//{
			//	in_label[*iter].insert(src, cur_distance);//cur_distance or + 1?
			//}
		}
		else if(min_order[cur_distance] == map_node_to_total_order[src]){
			if (!out_label[src].intersection(in_label[*iter], k))//if doesnot exist common node
			{
				//out_label[src].insert(*iter, cur_distance);
				in_label[*iter].insert(src, k-cur_distance);
			}
		}
	
		cur_distance++;
		cur_iters[cur_distance] = adjacency_list[ *cur_iters[cur_distance-1]].begin();
		
		cur_node = *cur_iters[cur_distance-1];
		c_path_set.insert(cur_node);
		c_path.push(cur_node);
		//cur_iters[cur_distance]++;

	}
}



//void construct_k_total_index(NODE_TYPE node_num, vector<map_distance_node_pair> &in_label, vector<map_distance_node_pair> &out_label, vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k)
void construct_k_total_index_dfs(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, NODE_TYPE src, vector<NODE_TYPE>& c_path, int k, NODE_TYPE cur_node, int cur_distance, set<NODE_TYPE>& c_path_set, vector<map_distance_node_pair> &in_label, vector<map_distance_node_pair> &out_label, NODE_TYPE mindis_in_curpath)// distance means length of path
{//only travel out-neibors, to get full label, should run adjacency_list reverse and reverse in out labels to get 
	if (cur_distance >= k)//
	{
		return;
	}
	cur_distance++;
	c_path.push_back(cur_node);
	//set<NODE_TYPE> new_c_path_set(c_path_set);
	c_path_set.insert(cur_node);
	for (auto iter = adjacency_list_double[cur_node].begin(); iter != adjacency_list_double[cur_node].end(); iter++)
	{//double directed graph, we can set adjacency_list_double to empty if we do not use it
		if (*iter < src || c_path_set.find(*iter) != c_path_set.end())
		{// there is a repeat node  or rank is smaller than src node
			continue;
		}
		int temp_min = mindis_in_curpath;
		if (*iter < mindis_in_curpath || mindis_in_curpath == -1)
		{//insert the new labels
			if (!out_label[src].intersection(in_label[*iter],k))//if doesnot exist common node
			{
				in_label[*iter].insert(src,cur_distance);
			}
			out_label[src].insert(*iter, cur_distance);
			temp_min = *iter;
		}
		construct_k_total_index_dfs(adjacency_list_double, adjacency_list, src, c_path, k, *iter, cur_distance, c_path_set, in_label, out_label, temp_min);
			//c_path_set.erase(*iter);
	}

	for (auto iter = adjacency_list[cur_node].begin(); iter != adjacency_list[cur_node].end(); iter++)
	{//double directed graph, we can set adjacency_list_double to empty if we do not use it
		if (*iter < src || c_path_set.find(*iter) != c_path_set.end())
		{// there is a repeat node  or rank is smaller than src node
			continue;
		}
		int temp_min = mindis_in_curpath;
		if (*iter < mindis_in_curpath || mindis_in_curpath == -1)
		{//insert the new labels

			if (!out_label[src].intersection(in_label[*iter], k))//if doesnot exist common node
			{
				in_label[*iter].insert(src, cur_distance);//cur_distance or + 1?
			}
			out_label[src].insert(*iter, cur_distance);
			temp_min = *iter;

		}
		construct_k_total_index_dfs(adjacency_list_double, adjacency_list, src, c_path, k, *iter, cur_distance, c_path_set, in_label, out_label, temp_min);

		//c_path_set.erase(*iter);
	}
	auto temp_del = c_path_set.find(cur_node);
	if (temp_del != c_path_set.end())
	{
		c_path_set.erase(temp_del);
	}
	c_path.pop_back();
}

set<NODE_TYPE> k_reach_from_node_s(NODE_TYPE s, NODE_TYPE node_num, vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k)
{//return k-reach nodes from source node s
	set<NODE_TYPE> result;
	set<NODE_TYPE> cur_propagate;
	int cur_dis = 0;
	result.insert(s);
	cur_propagate.insert(s);
	for (cur_dis = 0; cur_dis < k; cur_dis++)// since we add node in this iteration, so use < k
	{
		cur_dis++;
		set<NODE_TYPE> temp_proprogate;
		for (auto iter1 = cur_propagate.begin(); iter1 != cur_propagate.end(); iter1++)
		{
			for (auto iter2 = adjacency_list[*iter1].begin(); iter2 != adjacency_list[*iter1].end(); iter2++)
			{
				if (result.find(*iter2) != result.end())// if *iter2 is not in result set
				{
					result.insert(*iter2);
					temp_proprogate.insert(*iter2);
				}
			}
		}
		cur_propagate = temp_proprogate;
	}
	return result;
}


//extend algorithms to total order one

void dfs_with_total_order(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, map<NODE_TYPE, paths> &result, NODE_TYPE cur_node, vector<NODE_TYPE> c_path, set<NODE_TYPE> stop_nodes, int distance, int cur_distance, int& min_stop_disatance, set<NODE_TYPE> c_path_set, vector<map_distance_node_pair> &in_label, vector<map_distance_node_pair> &out_label, NODE_TYPE dst, set<NODE_TYPE> &label_nodes)// distance means length of path
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
	{//double directed graph, we can set adjacency_list_double to empty if we do not use it
		if (c_path_set.find(*iter) != c_path_set.end() || !(out_label[*iter].intersection(in_label[dst], distance)) )
		{// there is a repeat node or *iter cannot reach dst
			continue;
		}
		dfs_with_total_order(adjacency_list_double, adjacency_list, result, *iter, c_path, stop_nodes, distance, cur_distance + 1, min_stop_disatance, new_c_path_set,in_label, out_label,dst,label_nodes);
		//c_path_set.erase(*iter);
	}

	for (auto iter = adjacency_list[cur_node].begin(); iter != adjacency_list[cur_node].end(); iter++)
	{
		if (c_path_set.find(*iter) != c_path_set.end() )// there is a repeat node or *iter cannot reach dst
		{
			continue;
		}
		else if (label_nodes.find(*iter) != label_nodes.end() || label_nodes.find(dst) != label_nodes.end() ) {//*iter or dst are in label_nodes
			if (!(out_label[*iter].intersection(in_label[dst], distance)))// src cannot reach dst
			{
				continue;
			}
		}
		dfs_with_total_order(adjacency_list_double, adjacency_list, result, *iter, c_path, stop_nodes, distance, cur_distance + 1, min_stop_disatance, new_c_path_set,in_label,out_label,dst,label_nodes);
		//c_path_set.erase(*iter);
	}
}


void get_in_out_degrees(vector<int>& in_degrees, vector<int>& out_degrees, vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, int node_num)
{
	in_degrees.resize(node_num + 1);
	out_degrees.resize(node_num + 1);
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		for (auto iter = adjacency_list_double[i].begin(); iter != adjacency_list_double[i].end(); iter++)
		{
			in_degrees[*iter]++;
			out_degrees[*iter]++;
			in_degrees[i]++;
			out_degrees[i]++;
		}
		for (auto iter2 = adjacency_list[i].begin(); iter2 != adjacency_list[i].end(); iter2++)
		{
			in_degrees[*iter2]++;
			out_degrees[i]++;
		}

	}
}

//dfs hot point with total order
paths find_all_paths_with_hotpoints_dfs_with_total_order(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points, path_index &index, vector<map_distance_node_pair> &in_label, vector<map_distance_node_pair> &out_label, set<NODE_TYPE>& label_nodes)
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
	dfs_with_total_order(adjacency_list_double, adjacency_list, left_dfs, query_node1, c_path, stop_points, k, 0, min_stop_distance, c_path_set,in_label,out_label,query_node2,label_nodes);
	paths paths_without_hot;
	auto iter_left = left_dfs.find(query_node2);
	if (iter_left != left_dfs.end())
	{
		paths_without_hot.add_paths(iter_left->second);
		//paths_without_hot.output();
	}
	//paths_without_hot.output();
	set<NODE_TYPE> right_stop_points(hot_points);
	dfs_with_total_order(adjacency_list_double, adjacency_list_reverse, right_dfs, query_node2, c_path, right_stop_points, k - min_stop_distance, 0, right_stop_distance, c_path_set,in_label,out_label,query_node2,label_nodes);

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
			//left_hot_pahts_single_node.add_paths(temp_left);
			//right_hot_pahts_single_node.add_paths(temp_right);
			temp_left_distance = temp_left.get_min_distance();
			temp_right_distance = temp_right.get_min_distance();
			paths paths_bt_hot_points;
			cur_path c_path;
			index.find_paths_between_two_hot_nodes_index_without_cross_other_hotpoints(paths_bt_hot_points, c_path, left_node, right_node, (k - temp_left_distance - temp_right_distance), 0);
			//// cout << " paths bt hot points " << endl;
			//paths_bt_hot_points.output();
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
			index_new = temp_right.join(updated_edge, k);
			index_new = index_new.join(temp_left, k);
			index_new_paths.add_paths(index_new);

			//end of update	

		}
	}

	paths result;
	hot_paths.drop_path_length_less_than_k(1);// hot points path may contains path with length equal to 1
	result.add_paths(paths_without_hot);
	result.add_paths(hot_paths);
	//hot_paths.write_to_file("hint_step_by_step.txt");
	return result;

}

void consruct_I_in_out(NODE_TYPE  node_num, vector<map_distance_node_pair> &in_label, vector<map_distance_node_pair> &out_label, vector<map_distance_node_pair> &I_out, vector<map_distance_node_pair> &I_in)
{
	//construct I_in and I_out, put this into a seperate function
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		for (auto iter = in_label[i].v.begin(); iter != in_label[i].v.end(); iter++)
		{
			for (auto ele = iter->second->begin(); ele != iter->second->end(); ele++)
			{
				I_in[*ele].insert(i, iter->first);// insert i and distance information
			}
		}
		for (auto iter = out_label[i].v.begin(); iter != out_label[i].v.end(); iter++)
		{
			for (auto ele = iter->second->begin(); ele != iter->second->end(); ele++)
			{
				I_out[*ele].insert(i, iter->first);// insert i and distance information
			}
		}

	}
}


paths find_all_paths_with_hotpoints_dfs_with_total_order_dynamic(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  query_node1, NODE_TYPE  query_node2, NODE_TYPE  k, set<NODE_TYPE> hot_points, path_index &index, vector<map_distance_node_pair> &in_label, vector<map_distance_node_pair> &out_label, vector<map_distance_node_pair> &I_out, vector<map_distance_node_pair> &I_in, vector<NODE_TYPE>& map_node_to_total_order, int& cur_low_total_order, set<NODE_TYPE>& label_nodes)
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
	dfs_with_total_order(adjacency_list_double, adjacency_list, left_dfs, query_node1, c_path, stop_points, k, 0, min_stop_distance, c_path_set, in_label, out_label, query_node2,label_nodes);
	paths paths_without_hot;
	auto iter_left = left_dfs.find(query_node2);
	if (iter_left != left_dfs.end())
	{
		paths_without_hot.add_paths(iter_left->second);
		//paths_without_hot.output();
	}
	//paths_without_hot.output();
	set<NODE_TYPE> right_stop_points(hot_points);
	dfs_with_total_order(adjacency_list_double, adjacency_list_reverse, right_dfs, query_node2, c_path, right_stop_points, k - min_stop_distance, 0, right_stop_distance, c_path_set, in_label, out_label, query_node2,label_nodes);

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
			//left_hot_pahts_single_node.add_paths(temp_left);
			//right_hot_pahts_single_node.add_paths(temp_right);
			temp_left_distance = temp_left.get_min_distance();
			temp_right_distance = temp_right.get_min_distance();
			paths paths_bt_hot_points;
			cur_path c_path;
			index.find_paths_between_two_hot_nodes_index_without_cross_other_hotpoints(paths_bt_hot_points, c_path, left_node, right_node, (k - temp_left_distance - temp_right_distance), 0);
			//// cout << " paths bt hot points " << endl;
			//paths_bt_hot_points.output();
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
			//paths index_new;
			//vector<vector<NODE_TYPE> > edge;
			//vector<NODE_TYPE> query_edge;
			//query_edge.push_back(query_node2);
			//query_edge.push_back(query_node1);// the new updated edge is from query_node2 to query_node1(the direction)
			//edge.push_back(query_edge);
			//paths updated_edge(edge);
			//index_new = temp_right.join(updated_edge, k);
			//index_new = index_new.join(temp_left, k);
			//index_new_paths.add_paths(index_new);

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
		//index.push_back(index_new_paths);//already updated in hp index algorithm
		//adjacency_list[query_node2].push_back(query_node1);//if the query is from node1 to node2, then the new edge is node2->node1 to from a circle
		//adjacency_list_reverse[query_node1].push_back(query_node2);

		//update index
	
		//dynamic for k reach total order
		//new edge is node2->node1

		bool node2inlabel = false;
		bool node1inlabel = false;
		if (label_nodes.find(query_node2) != label_nodes.end())
		{
			node2inlabel = true;
		}
		if (label_nodes.find(query_node1) != label_nodes.end())
		{
			node1inlabel = true;
		}
		if (map_node_to_total_order[query_node1] < map_node_to_total_order[query_node2])
		{//for v in I_out[query_node2], we update labels
			if (node2inlabel && node1inlabel)
			{
				for (auto iter = I_out[query_node2].v.begin(); iter != I_out[query_node2].v.end(); iter++)
				{
					if (iter->first >= k)// skip k
					{
						continue;
					}
					set<NODE_TYPE> temp_set(*iter->second);
					for (auto ele = temp_set.begin(); ele != temp_set.end(); ele++)
					{
						NODE_TYPE t_node = *ele;
						//out_label[t_node].remove(*ele, iter->first);
						out_label[t_node].insert(query_node1, iter->first + 1);
						//I_out[query_node2].remove(*ele, iter->first);

						I_out[query_node1].insert(t_node, iter->first + 1);

					}
				}
			}
			//for querynode2 insert (1,1) and out_label[1]'s label with distance smaller than k
			if (node1inlabel) 
			{
				out_label[query_node2].insert(query_node1, 1);
			}
			for (auto iter = out_label[query_node1].v.begin(); iter != out_label[query_node1].v.end(); iter++)
			{
				if (iter->first >= k)
				{
					continue;
				}
				for (auto ele = iter->second->begin(); ele != iter->second->end(); ele++)
				{
					out_label[query_node2].insert(*ele, iter->first + 1);
				}
			}
			//for query_node1 insert (2,1) and in_label[2]'s label with distance smaller than k
			//if (node2inlabel)
			//{
			//	in_label[query_node1].insert(query_node2, 1);//no need since node2's total value is larger than node1
			//}
			for (auto iter = in_label[query_node2].v.begin(); iter != in_label[query_node2].v.end(); iter++)
			{
				if (iter->first >= k)
				{
					continue;
				}
				for (auto ele = iter->second->begin(); ele != iter->second->end(); ele++)
				{
					if (map_node_to_total_order[*ele] < map_node_to_total_order[query_node1] && label_nodes.find(*ele) != label_nodes.end())
					{
						in_label[query_node1].insert(*ele, iter->first + 1);
					}
				}
			}
		}
		else
		{//for v in I_in[query_node1], we update labels
			if (node2inlabel && node1inlabel) 
			{
				for (auto iter = I_in[query_node1].v.begin(); iter != I_in[query_node1].v.end(); iter++)
				{
					if (iter->first >= k)// skip k
					{
						continue;
					}
					set<NODE_TYPE> temp_set(*iter->second);
					for (auto ele = temp_set.begin(); ele != temp_set.end(); ele++)
					{
						NODE_TYPE t_node = *ele;
						//in_label[t_node].remove(*ele, iter->first);
						in_label[t_node].insert(query_node2, iter->first + 1);
						//I_in[query_node1].remove(*ele, iter->first);
						I_in[query_node2].insert(t_node, iter->first + 1);
						//update in_label

					}
				}
			}
			//for querynode1 insert (2,1) and in_label[2]'s label with distance smaller than k
			if (node2inlabel) 
			{
				in_label[query_node1].insert(query_node2, 1);
			}
			for (auto iter = in_label[query_node2].v.begin(); iter != in_label[query_node2].v.end(); iter++)
			{
				if (iter->first >= k)
				{
					continue;
				}
				for (auto ele = iter->second->begin(); ele != iter->second->end(); ele++)
				{
					in_label[query_node1].insert(*ele, iter->first + 1);
				}
			}
			//for query_node2 insert (1,1) and out_label[1]'s label with distance smaller than k
			//out_label[query_node2].insert(query_node1, 1);
			for (auto iter = out_label[query_node1].v.begin(); iter != out_label[query_node1].v.end(); iter++)
			{
				if (iter->first >= k)
				{
					continue;
				}
				for (auto ele = iter->second->begin(); ele != iter->second->end(); ele++)
				{
					if (map_node_to_total_order[*ele] < map_node_to_total_order[query_node2] && label_nodes.find(*ele) != label_nodes.end())
					{
						out_label[query_node2].insert(*ele, iter->first + 1);
					}
				}
			}
		}
	}
		//end of update

	paths result;
	hot_paths.drop_path_length_less_than_k(1);// hot points path may contains path with length equal to 1
	result.add_paths(paths_without_hot);
	result.add_paths(hot_paths);
	
	return result;

}

set<NODE_TYPE> find_hot_points_total_order(vector<NODE_TYPE >* adjacency_list_double, vector<NODE_TYPE >* adjacency_list, vector<NODE_TYPE >* adjacency_list_reverse, NODE_TYPE  node_num, NODE_TYPE  k, vector<NODE_TYPE>& map_node_to_total_order,int& cur_low_total_order, double threhold)
{
	vector<hot_degree> hot_points;
	//double threhold = 0.05;//this means we find top 50% nodes in degree rank as hot-points(this part can change to k-core value, truess value or centrality)
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		NODE_TYPE temp_degree = adjacency_list[i].size() + adjacency_list_double[i].size();
		if (temp_degree != 0)
		{
			hot_degree hd1(i, temp_degree);
			hot_points.push_back(hd1);
		}
	}
	stable_sort(hot_points.begin(), hot_points.end(), less<hot_degree>());
	NODE_TYPE end = (NODE_TYPE)(hot_points.size() * threhold);
	// cout << " total nodes is " << hot_points.size() << endl;
	// cout << " we choose " << end << endl;
	set<NODE_TYPE> result;
	for (NODE_TYPE i = 0; i < end; i++)
	{
		result.insert(hot_points[i].node);
		//change hot_points[i].node's map value to a new one
		NODE_TYPE origin_node = hot_points[i].node;
		int temp_value = map_node_to_total_order[origin_node];
		map_node_to_total_order[origin_node] = map_node_to_total_order[cur_low_total_order];
		map_node_to_total_order[cur_low_total_order] = temp_value;
		cur_low_total_order++;
	}
	// cout << "hot points : " << result.size() << endl;
	for (auto iter = result.begin(); iter != result.end(); iter++)
	{
		// cout << *iter << " ";
	}
	// cout << endl;
	return result;
}


//void static_k_reach_total_order_test()
//{
//	set<NODE_TYPE> label_nodes;
//	vector<NODE_TYPE> map_node_to_total_order;
//	int cur_low_total_order = 0;
//	int k = 3;
//	NODE_TYPE  node_num = 548552; //1007518272;//548552;//1234944765;
//								  //NODE_TYPE node_num = 91306;
//	vector<int> in_degrees, out_degrees;// we need to get in and out degrees
//	vector<NODE_TYPE >* adjacency_list = new vector<NODE_TYPE >[node_num + 1];
//	vector<NODE_TYPE >* adjacency_list_reverse = new vector<NODE_TYPE >[node_num + 1];
//	vector<NODE_TYPE >* adjacency_list_double = new vector<NODE_TYPE>[node_num + 1];//record the double directed adjacency for the undirected graph load 
//	load_graph_data("Amazon.txt", adjacency_list, adjacency_list_reverse);
//	vector<map_distance_node_pair> in_label;
//	vector<map_distance_node_pair> out_label;
//	in_label.resize(node_num);
//	out_label.resize(node_num);
//	set<NODE_TYPE> hot_points;
//	clock_t hpindex_start, hpindex_end;
//	hpindex_start = clock();
//	path_index index = construct_hot_point_index(adjacency_list, adjacency_list_reverse, node_num, k, hot_points);
//	hpindex_end = clock();
//	// cout << "end of hot point index" << double(hpindex_end - hpindex_start) / CLOCKS_PER_SEC << endl;
//	
//
//	get_in_out_degrees(in_degrees, out_degrees, adjacency_list_double, adjacency_list, node_num);
//
//	clock_t tindex_start, tindex_end;
//	tindex_start = clock();
//	construct_k_total_index(node_num, in_label, out_label, adjacency_list, adjacency_list_reverse,adjacency_list_double,in_degrees,out_degrees,k, map_node_to_total_order,cur_low_total_order, label_nodes);
//	tindex_end = clock();
//	// cout << "end of total order index" << double(tindex_end - tindex_start) / CLOCKS_PER_SEC << endl;
//	NODE_TYPE query_node1 = 2;
//	NODE_TYPE query_node2 = 411179;
//	paths result1;
//	paths result2;
//	clock_t start, end, end1;
//	start = clock();
//	result1 = find_all_paths_with_hotpoints_dfs(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k, hot_points, index);
//	end = clock();
//	result2 = find_all_paths_with_hotpoints_dfs_with_total_order(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k, hot_points, index, in_label, out_label,label_nodes);
//	end1 = clock();
//	// cout << "hot points dfs lasts for " << double(end - start) / CLOCKS_PER_SEC <<endl;
//	// cout << "total order labels dfs lasts for " << double(end1 - end) / CLOCKS_PER_SEC << endl;
//	result1.write_to_file("result1");
//	result2.write_to_file("result2");
//}

void dynamic_k_reach_total_order_test(int k, int option, int threshold1, double threshold2, long query_num, char* dataset, bool output_result, char* algorithm, bool random_query)
{//algorithms includes 
	char temp_algorithms[50][50] = { "induced_dag_meetnodes","dfsbaseline","double_bfs_baseline","induced_dag_dfs" };
	vector<NODE_TYPE> map_node_to_total_order;
	//int option = 0;// 0 means use dynamic files to test
	//double threshold1 = 0.0001;//hot points threshold
	//double threshold2 = 0.0010;//label nodes threshold
	//int k = 5;
							   // 1 means test with given query_node1 and query_node2
	int testi = 10;
	ofstream test;
	//test.open("Amazon_result\\hp_result" + to_string(testi));
	//test << "test " << endl;

	
	ifstream count_file;
	ifstream inEdges;
	NODE_TYPE  node_num;// = 548552; //121790;// 548552; //1007518272;//548552;//1234944765;
								  //NODE_TYPE node_num = 91306;
	string count_name = dataset;
	count_name += ".static";
	count_file.open(count_name);
	char tmp[BUFFER_LENTH];

	if (!count_file.is_open())
	{
		cout << "Error opening file : " << dataset << endl;
	}
	//configure the values of n and m 
	node_num = 0;
	int i = 0;
	NODE_TYPE x, y;
	while (!count_file.eof())
	{
		NODE_TYPE x, y;
		count_file.getline(tmp, BUFFER_LENTH);
		extractEdges(tmp, x, y);
		if (x > node_num)
		{
			node_num = x;
		}
		if (y > node_num)
		{
			node_num = y;
		}
		//adjacency_list[x].push_back(y);//directed graph, so we only push edge x y NODE_TYPE o x's adjacency list
		//adjacency_list_reverse[y].push_back(x);
		//adjacency_list[y].push_back(x);
	}
	//de_num++;
	node_num++;
	cout << node_num << " is the load node num " << endl;
	count_file.close();

	//if (strcmp(dataset, "Amazon.txt") == 0)
	//{
	//	node_num = 548552;
	//	//EDGE_COUNT = 2284995;
	//}
	//else if (strcmp(dataset, "gplus_combined.txtreoder") == 0)
	//{
	//	node_num = 121790;
	//}
	////else if (strcmp(dataset, "cit-Patents.txtreoder") == 0)
	////{
	////	node_num = 3774770;
	////}
	//else if (strcmp(dataset, "Slashdot0902.txtreoder") == 0)
	//{
	//	node_num = 82170;
	//}
	//else if (strcmp(dataset, "soc-Epinions1.txtreoder") == 0)
	//{
	//	node_num = 75880;
	//}
	//else if (strcmp(dataset, "soc-LiveJournal1.txtreoder") == 0)
	//{
	//	node_num = 4847580;
	//}
	//else if (strcmp(dataset, "soc-pokec-relationships.txtreoder") == 0)
	//{
	//	node_num = 1632810;
	//}
	//else if (strcmp(dataset, "WikiTalk.txtreoder") == 0)
	//{
	//	node_num = 2394390;
	//}
	//else if (strcmp(dataset, "twitter_social") == 0)
	//{
	//	node_num = 465020;
	//}
	//else if (strcmp(dataset, "com-friendster.ungraph.txt") == 0)
	//{
	//	node_num = 66608366;
	//}
	//else if (strcmp(dataset, "roadNet-CA.txt") == 0)
	//{
	//	node_num = 2966206;
	//}
	//else if (strcmp(dataset, "roadNet-PA.txt") == 0)
	//{
	//	node_num = 2096206;
	//}
	//else if (strcmp(dataset, "roadNet-TX.txt") == 0)
	//{
	//	node_num = 2466206;
	//}
	//else {
	//	// cout << "not match any datasets" << endl;
	//	node_num = 6608366;
	//	//return;
	//}
	//
	//dataset = "./dataset/" + dataset;

	vector<int> in_degrees, out_degrees;// we need to get in and out degrees
	vector<NODE_TYPE >* adjacency_list = new vector<NODE_TYPE >[node_num + 1];
	vector<NODE_TYPE >* adjacency_list_reverse = new vector<NODE_TYPE >[node_num + 1];
	vector<NODE_TYPE >* adjacency_list_double = new vector<NODE_TYPE>[node_num + 1];//record the double directed adjacency for the undirected graph load 





	//load_graph_data("gplus_combined.txtreoder.static", adjacency_list, adjacency_list_reverse);
	//inEdges.open("gplus_combined.txtreoder.dynamic");
	
	
	string temp_static(dataset);
	string temp_dynamic(dataset);

	if (!random_query)
	{
		temp_static += ".static";
		//string temp_static = "test_case.txt";

		temp_dynamic += ".dynamic";
		load_graph_data(temp_static, adjacency_list, adjacency_list_reverse);
	}
	else {
		temp_static += ".static";

		temp_dynamic += "_" + to_string(k) + "_randomquery";
		load_graph_data(temp_static, adjacency_list, adjacency_list_reverse);

	}
	//
	unordered_map<NODE_TYPE, set<NODE_TYPE>> adjacency_map;
	unordered_map<NODE_TYPE, set<NODE_TYPE>> adjacency_map_reverse;
	for(int i = 0; i <= node_num;i++)
	{
		bool first_time = true;
		for(auto iter_adj = adjacency_list[i].begin(); iter_adj != adjacency_list[i].end();iter_adj++)
		{
			if(first_time)
			{
				set<NODE_TYPE> temp_s;
				temp_s.insert((*iter_adj));
				adjacency_map.insert(std::make_pair(i, temp_s));
				first_time = false;
			}
			else
			{
				adjacency_map.find(i)->second.insert(*iter_adj);
			}
		}
	}

	for (int i = 0; i <= node_num; i++)
	{
		bool first_time = true;
		for (auto iter_adj = adjacency_list_reverse[i].begin(); iter_adj != adjacency_list_reverse[i].end(); iter_adj++)
		{
			if (first_time)
			{
				set<NODE_TYPE> temp_s;
				temp_s.insert((*iter_adj));
				adjacency_map_reverse.insert(std::make_pair(i, temp_s));
				first_time = false;
			}
			else
			{
				adjacency_map_reverse.find(i)->second.insert(*iter_adj);
			}
		}
	}

	inEdges.open(temp_dynamic);

	
	if (!inEdges.is_open())
	{
		cout << "Error opening file : " << temp_dynamic << endl;
	}
	vector<map_distance_node_pair> in_label;
	vector<map_distance_node_pair> out_label;
	vector<map_distance_node_pair> I_in;
	vector<map_distance_node_pair> I_out;
	int cur_low_total_order = 1;
	map_node_to_total_order.resize(node_num);
	for (int i = 0; i < node_num; i++)
	{
		map_node_to_total_order[i] = i;
	}
	in_label.resize(node_num);
	out_label.resize(node_num);
	I_in.resize(node_num);
	I_out.resize(node_num);
	//set<NODE_TYPE> hot_points = find_hot_points(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, k, threshold1);
	set<NODE_TYPE> hot_points = find_hot_points_top_t(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, k, threshold1);

	NODE_TYPE temp;
	//set<NODE_TYPE> hot_points;
	//ifstream hpf("test_hotpoints.txt");
	//while (!hpf.eof())
	//{
	//	hpf >> temp;
	//	// cout << temp << " in file" << endl;
	//	hot_points.insert(temp);
	//}
	
	set<NODE_TYPE> label_nodes = find_hot_points_total_order(adjacency_list,adjacency_list_double, adjacency_list_reverse, node_num, k,map_node_to_total_order,cur_low_total_order, threshold2);
	clock_t hpindex_start, hpindex_end;
	set<NODE_TYPE> empty_hot_points;

	hpindex_start = clock();
	//path_index index = construct_hot_point_index(adjacency_list, adjacency_list_reverse, node_num, k, test_hot_points);
	// cout << "start consturc hp index" << endl;
	//path_index index = construct_hot_point_index_dfs(adjacency_list, adjacency_list_reverse, k, hot_points);
	
	path_index index = construct_hot_point_index_dfs_using_new_algorithm(adjacency_map, k, hot_points,node_num);


	//paths temp_index =	index.find_paths_between_two_hot_nodes(501444, 502784, k);
	hpindex_end = clock();
	// cout << "end of hot point index" << double(hpindex_end - hpindex_start) / CLOCKS_PER_SEC << endl;


	get_in_out_degrees(in_degrees, out_degrees, adjacency_list_double, adjacency_list, node_num);

	clock_t tindex_start, tindex_end;
	tindex_start = clock();
	//construct_k_total_index(node_num, in_label, out_label, adjacency_list, adjacency_list_reverse, adjacency_list_double, in_degrees, out_degrees, k,map_node_to_total_order,cur_low_total_order,label_nodes);
	// cout << "end of total order index" << endl;
	consruct_I_in_out(node_num, in_label, out_label, I_out, I_in);
	//if (out_label[531483].intersection(in_label[2189], k))
	//{
	//	// cout << "true";
	//}
	tindex_end = clock();
	// cout << "end of total order index and I_in I_out" << double(tindex_end - tindex_start) / CLOCKS_PER_SEC << endl;


	clock_t spdindex_start, spdindex_end;
	spdindex_start = clock();
	spd_index spdindex;
	//spdindex.insert_hot_points_set(hot_points);
	//spdindex.construct_spd_index(adjacency_list, adjacency_list_reverse, node_num, k);
	spdindex_end = clock();


	NODE_TYPE query_node1 = 2;
	NODE_TYPE query_node2 = 6;
	paths result1;
	paths result2;
	paths result3[10];
	ofstream out_performance;
	ofstream repeat_rate;
	ofstream result_check;
	result_check.open("lines_result",ios::app);
	long long pruned_lines = 0;
	long long baseline_lines = 0;


	out_performance.open("runtime_"+to_string(k)+"_"+dataset+"_"+algorithm+"_"+to_string(threshold1));
	repeat_rate.open("repeat_rate");
	out_performance << "hot point(triple) index consturct time : " << double(hpindex_end - hpindex_start) / CLOCKS_PER_SEC << endl;
	out_performance << "total order(meet points) index consturct time : " << double(tindex_end - tindex_start) / CLOCKS_PER_SEC << endl;
	out_performance << "spd index consturct time : " << double(spdindex_end - spdindex_start) / CLOCKS_PER_SEC << endl;

	char buffer[255];
	//NODE_TYPE x, y;
	int edgeIndex = 0;
	clock_t start, end, end1, end2,start2, time_bf_pruned_path, time_bf_pruned_path_result5, t_construct_induced;
	set<NODE_TYPE> c_path_set2;
	vector<NODE_TYPE> c_path2;
	paths result5;
	pruned_subgraph_unordered_map_double_direction_only_meet_nodes_half_k temp_small_2;
	pruned_subgraph_unordered_map_double_direction temp_pruned_test;

	switch(option)
	{
	case 0:
		//get the dynamic edges and test for runtime and correctness
		
		out_performance << "edgeid \t " << algorithm << " query time " << endl;
		while (!inEdges.eof())
		{
			edgeIndex++;
			if (edgeIndex >= query_num)
			{
				break;
			}

			inEdges.getline(buffer, BUFFER_LENTH);
			extractEdges(buffer, query_node1, query_node2);
			if (edgeIndex <= 0)
			{
				continue;
			}
			cout << edgeIndex << "\t" << query_node1 << "\t" << query_node2 << "\t";
			int spd = bfs_return_spd(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			cout << spd << endl;
			set<NODE_TYPE> c_path_set;
			vector<NODE_TYPE> c_path;
			//result3[0]
			paths result4,result5;
			paths al_result;
			if (edgeIndex == 1)
			{
				al_result.clear_file("result" + to_string(k) + "_" + algorithm + "_" + dataset);
			}
			start = clock();
			t_construct_induced = clock();

			//choose algorithm

			if (strcmp(algorithm, "dfs") == 0)
			{
				dfs_find_k_paths(adjacency_list, query_node1, query_node1, query_node2, k, al_result, c_path_set, 0, c_path);
			}
			//find_path_within_length_constrain_only_reversed_shortestpathtree
			else if (strcmp(algorithm, "topk_dfs") == 0)
			{
				TopK_dfs temp_tes(adjacency_list, adjacency_list_reverse, query_node1, query_node2, k, node_num);
				temp_tes.construct_reverse_spd();
				al_result = temp_tes.find_path_within_length_constrain_only_reversed_shortestpathtree();
				//al_result.output();
			}
			else if (strcmp(algorithm, "topk_edbt") == 0)
			{
				TopK_dfs temp_tes(adjacency_list,adjacency_list_reverse,query_node1,query_node2,k,node_num);
				temp_tes.construct_reverse_spd();
				al_result = temp_tes.find_path_within_length_constrain();
				//al_result.output();
			}
			else if (strcmp(algorithm, "tdfs") == 0)
			{
				list_paths_theoryG(adjacency_list, adjacency_list_reverse, query_node1, query_node1, query_node2, k, al_result, c_path_set, 0, c_path);
			}
			else if (strcmp(algorithm, "tdfs_plus") == 0)
			{
				list_paths_theoryG_plus(adjacency_list, adjacency_list_reverse, query_node1, query_node1, query_node2, k, al_result, c_path_set, 0, c_path);
			}
			else if (strcmp(algorithm, "double_bfs") == 0)
			{
				al_result = double_direction_baseline_paths(adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k);
			}
			else if (strcmp(algorithm, "dfs_with_block") == 0)
			{
				CB_DFS dfs_block(adjacency_list,node_num);
				dfs_block.dfs_find_k_paths_with_block(query_node1,query_node1,query_node2,k,al_result,c_path_set,0,c_path);
				//al_result.drop_path_with_repeat_node();
			}
			else if (strcmp(algorithm, "dfs_with_block_induced") == 0)
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned2;
				temp_pruned2.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				CB_DFS_induced dfs_block(temp_pruned2.dag_min_induced_subgraph, node_num);
				//CB_DFS_induced dfs_block(temp_pruned2.dag_min_induced_subgraph, node_num,temp_pruned2.dst_distance);

				dfs_block.dfs_find_k_paths_with_block(query_node1, query_node1, query_node2, k, al_result, c_path_set, 0, c_path);
				//al_result.drop_path_with_repeat_node();
			}
			else if (strcmp(algorithm, "topk") == 0)
			{
				//al_result = ;
			}
			else if (strcmp(algorithm, "dfs_with_block_induced_spd_init") == 0)
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned2;
				temp_pruned2.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				//CB_DFS_induced dfs_block(temp_pruned2.dag_min_induced_subgraph, node_num);
				CB_DFS_induced dfs_block(temp_pruned2.dag_min_induced_subgraph, node_num, temp_pruned2.dst_distance);

				dfs_block.dfs_find_k_paths_with_block(query_node1, query_node1, query_node2, k, al_result, c_path_set, 0, c_path);
				//al_result.drop_path_with_repeat_node();
			}
			else if (strcmp(algorithm, "induced_dag_meetnodes") == 0)
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned3;
				temp_pruned3.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				temp_pruned3.find_all_meet_nodes_in_induced_subgraph(node_num, k, query_node1, query_node2);
				al_result = temp_pruned3.find_all_k_pahts_dfs_write_number_left_right_path(node_num, k, query_node1, query_node2,algorithm,dataset,edgeIndex);
			}
			else if (strcmp(algorithm, "induced_dag_dfs") == 0)
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned2;
				temp_pruned2.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				al_result = temp_pruned2.find_k_path_using_index(temp_pruned2.dag_min_induced_subgraph, node_num, k, query_node2, query_node1);
			}
			else if (strcmp(algorithm, "triple_join") == 0)
			{
				pruned_subgraph_unordered_map_double_direction_triple_join temp_triple;
				temp_triple.construct_pruned_subgraph_with_find_leftnodes_rightnodes_commonnodes(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				al_result = temp_triple.find_all_k_pahts_triple_join(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);

			}
			else if (strcmp(algorithm, "double_bfs_cut") == 0)
			{
				al_result.clear();
				pruned_subgraph_unordered_map_double_direction temp_pruned4;
				set<NODE_TYPE> stop_nodes;
				temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				unordered_map<NODE_TYPE, set<NODE_TYPE>> induced_subgraph_reverse = temp_pruned4.get_new_reverse_induced_subgraph(temp_pruned4.dag_min_induced_subgraph);
				short_cut_index short_index;
				short_index.construct_src_dst_short_paths(query_node1, query_node2, 3, temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse);

				/* test code for short cuts
				//paths temp = temp_pruned4.find_k_path_using_index_reverse(induced_subgraph_reverse, node_num, 2, query_node1, query_node2);
				//bool temp_bool = true;
				//paths temp2 = short_index.find_short_cuts_using_index(query_node1, query_node2, 2, temp_bool);
				*/
				al_result = temp_pruned4.find_k_paths_recursive_bfs_middle_nodes_only_bfs(temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse, node_num, query_node1, query_node2, k, stop_nodes, 1, short_index);

			}
			else if (strcmp(algorithm, "middle_nodes_cut_algirthm") == 0)//cut node order from spd small to large
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned4;

				set<NODE_TYPE> stop_nodes;
				short_cut_index short_index;
				temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				unordered_map<NODE_TYPE, set<NODE_TYPE>> induced_subgraph_reverse = temp_pruned4.get_new_reverse_induced_subgraph(temp_pruned4.dag_min_induced_subgraph);

				//short_index.construct_src_dst_short_paths(query_node1, query_node2, 3, temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse);
				al_result = temp_pruned4.find_k_paths_recursive_middle_cut_write_pathnum(temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse, node_num, query_node1, query_node2, k, stop_nodes, 1, short_index, algorithm, dataset, edgeIndex);

			}
			else if (strcmp(algorithm, "classical_nodes_cut_algirthm") == 0)//cut node order from spd small to large
			{
				//pruned_subgraph_unordered_map_double_direction temp_pruned4;

				//set<NODE_TYPE> stop_nodes;
				//short_cut_index short_index;
				//temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				//unordered_map<NODE_TYPE, set<NODE_TYPE>> induced_subgraph_reverse = temp_pruned4.get_new_reverse_induced_subgraph(temp_pruned4.dag_min_induced_subgraph);

				//al_result = temp_pruned4.find_k_paths_recursive_bfs_middle_nodes_write_pathnum(temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse, node_num, query_node1, query_node2, k, stop_nodes, 1, short_index, algorithm, dataset, edgeIndex);

			}
			else if (strcmp(algorithm, "double_bfs_cut_with_dfs") == 0)//cut node order from spd small to large
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned4;

				set<NODE_TYPE> stop_nodes;
				short_cut_index short_index;
				temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				unordered_map<NODE_TYPE, set<NODE_TYPE>> induced_subgraph_reverse = temp_pruned4.get_new_reverse_induced_subgraph(temp_pruned4.dag_min_induced_subgraph);

				//short_index.construct_src_dst_short_paths(query_node1, query_node2, 3, temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse);
				al_result = temp_pruned4.find_k_paths_recursive_bfs_middle_nodes_write_pathnum(temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse, node_num, query_node1, query_node2, k, stop_nodes, 1, short_index, algorithm, dataset, edgeIndex);

			}
			else if (strcmp(algorithm, "double_bfs_cut_with_dfs_random_order") == 0)//cut node order from spd small to large
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned4;

				set<NODE_TYPE> stop_nodes;
				short_cut_index short_index;
				temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				t_construct_induced = clock();
				ofstream runtime_construct;
				if (edgeIndex == 1)
				{
					runtime_construct.open("runtime_" + to_string(k) + "_" + dataset + "_" + algorithm + "_" + to_string(threshold1)+"_construct_induced");
				}
				else {
					runtime_construct.open("runtime_" + to_string(k) + "_" + dataset + "_" + algorithm + "_" + to_string(threshold1) + "_construct_induced", ios::app);

				}
				runtime_construct << edgeIndex << " \t " << double(t_construct_induced - start) / CLOCKS_PER_SEC << endl;
				unordered_map<NODE_TYPE, set<NODE_TYPE>> induced_subgraph_reverse = temp_pruned4.get_new_reverse_induced_subgraph(temp_pruned4.dag_min_induced_subgraph);

				//short_index.construct_src_dst_short_paths(query_node1, query_node2, 3, temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse);
				al_result = temp_pruned4.find_k_paths_recursive_bfs_middle_nodes_node_order_random(temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse, node_num, query_node1, query_node2, k, stop_nodes, 1, short_index);

			}
			else if (strcmp(algorithm, "double_bfs_cut_with_dfs_cost_model_order") == 0)//cut node order from spd small to large
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned4;

				set<NODE_TYPE> stop_nodes;
				short_cut_index short_index;
				temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				t_construct_induced = clock();
				ofstream runtime_construct;
				if (edgeIndex == 1)
				{
					runtime_construct.open("runtime_" + to_string(k) + "_" + dataset + "_" + algorithm + "_" + to_string(threshold1) + "_construct_induced");
				}
				else {
					runtime_construct.open("runtime_" + to_string(k) + "_" + dataset + "_" + algorithm + "_" + to_string(threshold1) + "_construct_induced", ios::app);

				}
				runtime_construct << edgeIndex << " \t " << double(t_construct_induced - start) / CLOCKS_PER_SEC << endl;
				unordered_map<NODE_TYPE, set<NODE_TYPE>> induced_subgraph_reverse = temp_pruned4.get_new_reverse_induced_subgraph(temp_pruned4.dag_min_induced_subgraph);

				//short_index.construct_src_dst_short_paths(query_node1, query_node2, 3, temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse);
				al_result = temp_pruned4.find_k_paths_recursive_bfs_middle_nodes_cut_node_order_using_costmodel(temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse, node_num, query_node1, query_node2, k, stop_nodes, 1, short_index);

			}
			else if (strcmp(algorithm, "double_bfs_cut_with_dfs_withoutdag") == 0)//cut node order from spd small to large
			{
				pruned_subgraph_unordered_map_double_direction temp_pruned4;
				set<NODE_TYPE> stop_nodes;
				short_cut_index short_index;
				//short_index.construct_src_dst_short_paths(query_node1, query_node2, 3, temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse);
				al_result = temp_pruned4.find_k_paths_recursive_bfs_middle_nodes_write_pathnum(adjacency_map, adjacency_map_reverse, node_num, query_node1, query_node2, k, stop_nodes, 1, short_index,algorithm,dataset,edgeIndex);
				
			}
			else if (strcmp(algorithm, "hp_index") == 0)//cut node order from spd small to large
			{
				double update_time = 0;
				al_result = find_all_paths_with_hotpoints_dfs_get_update_time_new_algor(adjacency_list_double, adjacency_list, adjacency_map, adjacency_map_reverse, node_num, query_node1, query_node2, k, hot_points, index, update_time);
				//need to delete update_time if want to get only query time
			}

			else {
				// cout << "not match any algorithms" << endl;
				node_num = 0;
				return;
			}

			end = clock();
			//al_result.drop_path_with_repeat_node();
			//al_result.sort_by_string_order();
			//al_result.drop_repeat_path();
			
			//adjacency_list[query_node2].push_back(query_node1);//if the query is from node1 to node2, then the new edge is node2->node1 to from a circle
			//adjacency_list_reverse[query_node1].push_back(query_node2);
			
			if (output_result)
			{
				al_result.write_to_file_append_edgeID("result" + to_string(k) + "_"+ algorithm  +"_" + dataset, edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
			}
			else
			{
				al_result.write_to_file_append_edgeID_and_result_size("result" + to_string(k) + "_" + algorithm + "_" + dataset, edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
			}
			out_performance << edgeIndex << " \t " << double(end - t_construct_induced) / CLOCKS_PER_SEC  << endl;
		}
		break;
	case 1:
		out_performance << "edgeid \t " << " pruned bfs query time \t" << "\hp index \t" << "total order \t" << "dfs baseline \t" << endl;
		while (!inEdges.eof())
		{
			edgeIndex++;
			if (edgeIndex >= query_num)
			{
				break;
			}
			// cout << edgeIndex << endl;
			inEdges.getline(buffer, BUFFER_LENTH);
			extractEdges(buffer, query_node1, query_node2);

			if (DEBUG_ALL_QUERYS) {
				ofstream temp;
				temp.open("allQuerys", ios::app);
				if (!temp)
				{
					// cout << "allQuerys can't open" << endl;
					abort();
				}
				temp << edgeIndex << " edge id " << endl;
				temp.close();
			}
			set<NODE_TYPE> c_path_set;
			vector<NODE_TYPE> c_path;
			//result3[0]
			paths result4, result5;
			start2 = clock();
			//result5 = find_k_paths_between_two_nodes_spd(adjacency_list_double, adjacency_list, adjacency_list_reverse, k, query_node1, query_node2,spdindex);
			//double_direction_baseline(adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k);


			//pruned_subgraph_unordered_map_double_direction temp_pruned;
			//temp_pruned.construct_pruned_subgraph(adjacency_list,adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//temp_pruned.join_left_right_index_into_left();
			//result5 = temp_pruned.find_k_path_using_index_reverse(temp_pruned.reverse_adjacency_in_subgraph_left,node_num,k,query_node1,query_node2);
			//
			
			//pruned_subgraph_unordered_map_double_direction_triple_join temp_triple;
			//temp_triple.construct_pruned_subgraph_with_find_leftnodes_rightnodes_commonnodes(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//result5 = temp_triple.find_all_k_pahts_triple_join(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//

			//pruned_subgraph_unordered_map_double_direction temp_pruned2;
			//temp_pruned2.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//result5 = temp_pruned2.find_k_path_using_index(temp_pruned2.dag_min_induced_subgraph, node_num, k, query_node2, query_node1);
			
			pruned_subgraph_unordered_map_double_direction temp_pruned3;
			temp_pruned3.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			temp_pruned3.find_all_meet_nodes_in_induced_subgraph(node_num, k, query_node1, query_node2);
			result5 = temp_pruned3.find_all_k_pahts_dfs(node_num, k, query_node1, query_node2);

			//without induced subgraph
			//temp_pruned3.find_all_meet_nodes_in_induced_subgraph_with_provided_induced_sbugraph(adjacency_map, adjacency_map_reverse, node_num, k, query_node1, query_node2);
			//result5 = temp_pruned3.find_all_k_pahts_dfs_with_provided_subgraph(adjacency_map, adjacency_map_reverse, node_num, k, query_node1, query_node2);

			





			time_bf_pruned_path_result5 = clock();
			//result5.drop_path_with_repeat_node();
			//result5.sort_by_string_order();
			//result5.drop_repeat_path();
			// cout << "end of meetnodes" << endl;
			//pruned_subgraph_unordered_map final_pruned;
			//final_pruned.construct_pruned_subgraph_by_pruned_subgraph(temp_pruned, node_num, k, query_node2, query_node1);
			//result5 = final_pruned.find_k_path_using_index(node_num, k, query_node2, query_node1);
			end2 = clock();
			//pruned_subgraph_unordered_map_double_direction temp_pruned2;
			//temp_pruned2.construct_pruned_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//temp_pruned2.join_right_left_into_right();
			//result4 = temp_pruned2.find_k_path_using_index(temp_pruned2.reverse_adjacency_in_subgraph_right, node_num, k, query_node2, query_node1);




			pruned_subgraph_unordered_map_double_direction temp_small;
			bool left = true;
			//random vertices order start
			/*pruned_subgraph_unordered_map_double_direction temp_pruned_random;
			set<NODE_TYPE> stop_nodes_random;
			short_cut_index short_index_random;

			result4 = temp_pruned_random.find_k_paths_recursive_bfs_middle_nodes_node_order_random(adjacency_map, adjacency_map_reverse, node_num, query_node1, query_node2, k, stop_nodes_random, 1, short_index_random);
			sort(result4.path.begin(), result4.path.end());
			*/

			//random vertices order end



			//temp_small.construct_pruned_subgraph_equal_double_degrees(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2,in_degrees, out_degrees,left);
			////if (temp_small.reverse_adjacency_in_subgraph_left.size() >= temp_small.reverse_adjacency_in_subgraph_right.size())
			//if(!left)
			//{
			//	temp_small.join_left_right_index_into_left();
			//	result4 = temp_small.find_k_path_using_index_reverse(temp_small.reverse_adjacency_in_subgraph_left, node_num, k, query_node1, query_node2);

			//}
			//else {
			//	temp_small.join_right_left_into_right();
			//	result4 = temp_small.find_k_path_using_index(temp_small.reverse_adjacency_in_subgraph_right, node_num, k, query_node2, query_node1);

			//}


			//dag minimum induced subgraph




			//dfs
			//set<NODE_TYPE> stop_points2;
			//map<NODE_TYPE, paths> left_dfs;
			//map<NODE_TYPE, paths> right_dfs;
			//DISTANCE_TYPE min_stop_distance;
			//vector<NODE_TYPE> c_path_1;
			//set<NODE_TYPE> c_path_set_1;
			//dfs_ex_node2(adjacency_list_double, adjacency_list, left_dfs, query_node1, c_path_1, stop_points2, k, 0, min_stop_distance, c_path_set_1, query_node2);
			////paths paths_without_hot;
			//auto iter_left = left_dfs.find(query_node2);
			//if (iter_left != left_dfs.end())
			//{
			//	result4.add_paths(iter_left->second);
			//	//paths_without_hot.output();
			//}

			dfs_find_k_paths(adjacency_list, query_node1, query_node1, query_node2, k, result4, c_path_set, 0, c_path);
			/*induced dfs
			pruned_subgraph_unordered_map_double_direction temp_pruned5;
			temp_pruned5.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			result4 = temp_pruned5.find_k_path_using_index(temp_pruned5.dag_min_induced_subgraph, node_num, k, query_node2, query_node1);
			*/
			


			//pruned_subgraph_unordered_map_double_direction temp_pruned3;
			//temp_pruned3.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//temp_pruned3.find_all_meet_nodes_in_induced_subgraph(node_num, k, query_node1, query_node2);
			//result4 = temp_pruned3.find_all_k_pahts_dfs(node_num, k, query_node1, query_node2);
			//result4.sort_by_string_order();



			//result4 = temp_pruned2.find_k_path_using_index_reverse(temp_pruned2.reverse_adjacency_in_subgraph_left, node_num, k, query_node1, query_node2);



			start = clock();			
			// cout << "end of dfs" << double(start - end2)/ CLOCKS_PER_SEC << endl;

			result1.clear();
			//pruned_subgraph_unordered_map_double_direction temp_pruned6;
			//temp_pruned6.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//temp_pruned6.join_left_right_index_into_left();
			//result1 = temp_pruned6.find_k_path_using_index_reverse(temp_pruned6.reverse_adjacency_in_subgraph_left, node_num, k, query_node1, query_node2);
			//pruned_subgraph_unordered_map_double_direction temp_pruned4;
			//set<NODE_TYPE> stop_nodes;
			//temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//unordered_map<NODE_TYPE, set<NODE_TYPE>> induced_subgraph_reverse = temp_pruned4.get_new_reverse_induced_subgraph(temp_pruned4.dag_min_induced_subgraph);
			//short_cut_index short_index;
			//short_index.construct_src_dst_short_paths(query_node1, query_node2, 3, temp_pruned4.dag_min_induced_subgraph, induced_subgraph_reverse);

			/* test code for short cuts
			//paths temp = temp_pruned4.find_k_path_using_index_reverse(induced_subgraph_reverse, node_num, 2, query_node1, query_node2);
			//bool temp_bool = true;
			//paths temp2 = short_index.find_short_cuts_using_index(query_node1, query_node2, 2, temp_bool);
			*/
			
			//result1 = temp_pruned4.find_k_paths_recursive_bfs_middle_nodes( temp_pruned4.dag_min_induced_subgraph,induced_subgraph_reverse, node_num, query_node1, query_node2, k, stop_nodes,1, short_index);
			
			//result1 = temp_pruned4.find_k_paths_recursive_bfs_middle_nodes(adjacency_map,adjacency_map_reverse, node_num, query_node1, query_node2, k, stop_nodes,1,short_index);
			
			//long before_remove = result1.path.size();
			
			//sort(result1.path.begin(), result1.path.end());
			/*auto it = unique(result1.path.begin(), result1.path.end());
			result1.path.erase(it, result1.path.end());*/
			//long after_remove = result1.path.size();
			//ofstream temp;
			//temp.open("repeat_ours", ios::app);
			//if (!temp)
			//{
			//	// cout << "repeat_ours can't open" << endl;
			//	abort();
			//}
			////temp << query_node1 << "\t" << query_node2 << " \t " << k << endl;
			//temp << before_remove << "\t" << after_remove << "\t" << double(before_remove)/after_remove << "\t for our al" << endl;
			//// cout << before_remove << "\t" << after_remove << "\t" << double(before_remove) / after_remove << "\t for our al" << endl;
			//temp.close();

			/*result1.sort_by_string_order();
			result1.drop_repeat_path();*/

			/*pruned_subgraph_unordered_map_double_direction temp_pruned4;
			temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//temp_pruned4.find_all_meet_nodes_in_induced_subgraph(node_num, k, query_node1, query_node2);
			temp_pruned4.find_all_meet_nodes_longspd_in_induced_subgraph(node_num, k, query_node1, query_node2);
			result1 = temp_pruned4.find_all_k_pahts_dfs_from_middle(node_num, k, query_node1, query_node2);
			*/
			
			
			//result1.drop_path_with_repeat_node();
			//// cout << "end of result1" << endl;
			//result1.sort_by_string_order();
			//result1.drop_repeat_path();
			//// cout << "end of result1 drop" << endl;


			//result1 = double_direction_baseline_paths(adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k);
			
			
			//pruned_subgraph_unordered_map temp_pruned2;
			//temp_pruned2.construct_pruned_subgraph(adjacency_list, node_num, k, query_node1, query_node2);
			//temp_pruned2.reverse();
			//result1 = temp_pruned2.find_k_path_using_index(node_num, k, query_node2, query_node1);
			/*pruned_subgraph_unordered_map_double_direction_triple_join temp_pruned2;
			temp_pruned2.construct_pruned_subgraph_with_find_leftnodes_rightnodes_commonnodes(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			result1 = temp_pruned2.find_all_k_pahts_triple_join(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);// temp_pruned2.find_k_path_using_index(temp_pruned2.reverse_adjacency_in_subgraph_right, node_num, k, query_node2, query_node1);
			result1.drop_path_with_repeat_node();
			result1.sort_by_string_order();
			result1.drop_repeat_path();
			*/
			double update_time = 0;

			//result1 = find_all_paths_with_hotpoints_dfs_get_update_time(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k, hot_points, index, update_time);
			result1 = find_all_paths_with_hotpoints_dfs_get_update_time_new_algor(adjacency_list_double, adjacency_list,adjacency_map, adjacency_map_reverse,node_num, query_node1, query_node2, k, hot_points, index, update_time);
			
			//find_all_paths_with_hotpoints_dfs_get_update_time
			//result1 = find_all_paths_with_hotpoints_dfs(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k, hot_points, index);
			end = clock();
			//bfs meetnodes 
			//pruned_subgraph_unordered_map_double_direction_only_meet_nodes_half_k temp_small_1;
			//temp_small_1.construct_pruned_subgraph_with_meetnodes(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			//result2 = temp_small_1.find_all_k_pahts_dfs(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			result2.clear();
			//dfs_find_k_paths(adjacency_list, query_node1, query_node1, query_node2, k, result2, c_path_set, 0, c_path);
			
			
			//result2.drop_path_with_repeat_node();
			//result2.sort_by_string_order();
			//result2.drop_repeat_path();


			//double_direction_baseline(adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k);
			//bool left = true;
			//temp_small_1.construct_pruned_subgraph_equal_double_degrees(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2, in_degrees, out_degrees, left);
			////if (temp_small.reverse_adjacency_in_subgraph_left.size() >= temp_small.reverse_adjacency_in_subgraph_right.size())
			//if (left)
			//{
			//	temp_small_1.join_left_right_index_into_left();
			//	result2 = temp_small_1.find_k_path_using_index_reverse(temp_small_1.reverse_adjacency_in_subgraph_left, node_num, k, query_node1, query_node2);

			//}
			//else {
			//	temp_small_1.join_right_left_into_right();
			//	result2 = temp_small_1.find_k_path_using_index(temp_small_1.reverse_adjacency_in_subgraph_right, node_num, k, query_node2, query_node1);

			//}
			time_bf_pruned_path = clock();
			int bf = result2.get_path().size();
			//result2.drop_path_with_repeat_node();
			int pruned_repeat_node = result2.get_path().size();




			result2.sort_by_string_order();
			result2.drop_repeat_path();
			//result2.sort_by_string_order();
			int after = result2.get_path().size();
			repeat_rate << "result2 : \t" << edgeIndex << " \t " << query_node2 << " : " << query_node1 << " : \t " << to_string(double(bf) / after) << " : " << bf << " : " << pruned_repeat_node << " : " << after << endl;
			//result2 = find_all_paths_with_hotpoints_dfs_with_total_order_dynamic(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k, empty_hot_points, index, in_label, out_label, I_out, I_in,map_node_to_total_order,cur_low_total_order,label_nodes);
			end1 = clock();

			result4.drop_path_with_repeat_node();
			result4.sort_by_string_order();
			result4.drop_repeat_path();

			result5.drop_path_with_repeat_node();
			result5.sort_by_string_order();
			result5.drop_repeat_path();

			//adjacency_list[query_node2].push_back(query_node1);//if the query is from node1 to node2, then the new edge is node2->node1 to from a circle
			//adjacency_list_reverse[query_node1].push_back(query_node2);
			result1.drop_path_with_repeat_node();
			result1.drop_path_not_start_from_nodeS(query_node1);
			//update adjacency_map and adjacency_map_reverse
			/*if(adjacency_map.find(query_node2) == adjacency_map.end())
			{
				set<NODE_TYPE> temp_adj;
				temp_adj.insert((query_node1));
			}
			else
			{
				adjacency_map.find(query_node2)->second.insert(query_node1);
			}
			if (adjacency_map_reverse.find(query_node1) == adjacency_map_reverse.end())
			{
				set<NODE_TYPE> temp_adj;
				temp_adj.insert((query_node2));
			}
			else
			{
				adjacency_map_reverse.find(query_node1)->second.insert(query_node2);
			}
			*/


			// cout << "spd dfs lasts for " << double(end2 - start2) / CLOCKS_PER_SEC << endl;
			// cout << "hot points dfs lasts for " << double(end - start) / CLOCKS_PER_SEC - update_time << endl;
			// cout << "total order labels dfs lasts for " << double(end1 - end) / CLOCKS_PER_SEC << endl;
			// cout << "baseline dfs lasts for " << double(start - end2) / CLOCKS_PER_SEC << endl;




			result2.sort_by_string_order();
			result4.sort_by_string_order();
			result1.sort_by_string_order();
			if (output_result)
			{
				result1.write_to_file_append_edgeID("hp_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
				result2.write_to_file_append_edgeID("total_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
				result4.write_to_file_append_edgeID("dfs_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);
				result5.write_to_file_append_edgeID("bfs_pruned_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);
			}
			else
			{
				result1.write_to_file_append_edgeID_and_result_size("hp_result" + to_string(k) + "_" + dataset+"_"+to_string(threshold1), edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
				result2.write_to_file_append_edgeID_and_result_size("total_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
				result4.write_to_file_append_edgeID_and_result_size("dfs_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);
				result5.write_to_file_append_edgeID_and_result_size("bfs_pruned_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);
			}
			pruned_lines += result5.path.size();
			baseline_lines += result1.path.size();
			out_performance << edgeIndex << " \t " << double(end2 - start2) / CLOCKS_PER_SEC << "\t" << double(end - start) / CLOCKS_PER_SEC - update_time << "\t " << double(end1 - end) / CLOCKS_PER_SEC << "\t" << double(start - end2) / CLOCKS_PER_SEC << "\t" << double(time_bf_pruned_path - end) / CLOCKS_PER_SEC << "\t" << double(time_bf_pruned_path_result5 - start2) / CLOCKS_PER_SEC << endl;
			//if (pruned_lines != baseline_lines)
			//{
			//	result_check << to_string(k) << "\t" << edgeIndex << "\t pruned_lines : " << pruned_lines << "\t baseline_lines : " << baseline_lines << endl;
			//	result1.write_to_file_append_edgeID("dfs_result_error", edgeIndex, query_node1, query_node2);
			//	result5.write_to_file_append_edgeID("bfs_pruned_result_error", edgeIndex, query_node1, query_node2);
			//	break;
			//}
			//break;
		}
		break;
	case 2:
		out_performance << "edgeid \t " << " pruned bfs query time \t" << "\hp index \t" << "total order \t" << "dfs baseline \t" << endl;
		while (!inEdges.eof())
		{
			edgeIndex++;
			if (edgeIndex >= query_num)
			{
				break;
			}
			// cout << edgeIndex << endl;
			random_query_pair(query_node1, query_node2, node_num);
			// cout << query_node1 << "\t" << query_node2 << endl;
			//inEdges.getline(buffer, BUFFER_LENTH);
			//extractEdges(buffer, query_node2, query_node1);
			//if (edgeIndex <= 28)
			//{
				//continue;
			//}
			if (DEBUG_ALL_QUERYS) {
				ofstream temp;
				temp.open("allQuerys", ios::app);
				if (!temp)
				{
					// cout << "allQuerys can't open" << endl;
					abort();
				}
					temp << edgeIndex << " edge id " << endl;
					temp.close();
				}
				set<NODE_TYPE> c_path_set;
				vector<NODE_TYPE> c_path;
				//result3[0]
				paths result4, result5;
				start2 = clock();

				pruned_subgraph_unordered_map_double_direction temp_pruned3;
				temp_pruned3.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				temp_pruned3.find_all_meet_nodes_in_induced_subgraph(node_num, k, query_node1, query_node2);
				result5  = temp_pruned3.find_k_path_using_index(temp_pruned3.dag_min_induced_subgraph, node_num, k, query_node2, query_node1);


				//without induced subgraph
				//temp_pruned3.find_all_meet_nodes_in_induced_subgraph_with_provided_induced_sbugraph(adjacency_map, adjacency_map_reverse, node_num, k, query_node1, query_node2);
				//result5 = temp_pruned3.find_all_k_pahts_dfs_with_provided_subgraph(adjacency_map, adjacency_map_reverse, node_num, k, query_node1, query_node2);







				time_bf_pruned_path_result5 = clock();
				end2 = clock();



				pruned_subgraph_unordered_map_double_direction temp_small;
				bool left = true;
				//random vertices order start
				/*pruned_subgraph_unordered_map_double_direction temp_pruned_random;
				set<NODE_TYPE> stop_nodes_random;
				short_cut_index short_index_random;

				result4 = temp_pruned_random.find_k_paths_recursive_bfs_middle_nodes_node_order_random(adjacency_map, adjacency_map_reverse, node_num, query_node1, query_node2, k, stop_nodes_random, 1, short_index_random);
				sort(result4.path.begin(), result4.path.end());
				*/

				//random vertices order end




				start = clock();
				// cout << "end of dfs" << double(start - end2) / CLOCKS_PER_SEC << endl;

				result1.clear();

				pruned_subgraph_unordered_map_double_direction temp_pruned4;
				temp_pruned4.construct_pruned_dag_min_subgraph(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
				//temp_pruned4.find_all_meet_nodes_in_induced_subgraph(node_num, k, query_node1, query_node2);
				temp_pruned4.find_all_meet_nodes_in_induced_subgraph(node_num, k, query_node1, query_node2);
				result1 = temp_pruned3.find_all_k_pahts_dfs_write_number_left_right_path(node_num, k, query_node1, query_node2, algorithm, dataset, edgeIndex);







				double update_time = 0;

				//result1 = find_all_paths_with_hotpoints_dfs_get_update_time(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k, hot_points, index, update_time);
				
				//hp index 
				//result1 = find_all_paths_with_hotpoints_dfs_get_update_time_new_algor(adjacency_list_double, adjacency_list, adjacency_map, adjacency_map_reverse, node_num, query_node1, query_node2, k, hot_points, index, update_time);

				end = clock();
				result2.clear();

				time_bf_pruned_path = clock();
				int bf = result2.get_path().size();
				//result2.drop_path_with_repeat_node();
				int pruned_repeat_node = result2.get_path().size();




				result2.sort_by_string_order();
				result2.drop_repeat_path();
				//result2.sort_by_string_order();
				int after = result2.get_path().size();
				repeat_rate << "result2 : \t" << edgeIndex << " \t " << query_node2 << " : " << query_node1 << " : \t " << to_string(double(bf) / after) << " : " << bf << " : " << pruned_repeat_node << " : " << after << endl;
				//result2 = find_all_paths_with_hotpoints_dfs_with_total_order_dynamic(adjacency_list_double, adjacency_list, adjacency_list_reverse, node_num, query_node1, query_node2, k, empty_hot_points, index, in_label, out_label, I_out, I_in,map_node_to_total_order,cur_low_total_order,label_nodes);
				end1 = clock();

				result4.drop_path_with_repeat_node();
				result4.sort_by_string_order();
				result4.drop_repeat_path();

				result5.drop_path_with_repeat_node();
				result5.sort_by_string_order();
				result5.drop_repeat_path();

				//adjacency_list[query_node2].push_back(query_node1);//if the query is from node1 to node2, then the new edge is node2->node1 to from a circle
				//adjacency_list_reverse[query_node1].push_back(query_node2);
				result1.drop_path_with_repeat_node();
				result1.drop_path_not_start_from_nodeS(query_node1);
				//update adjacency_map and adjacency_map_reverse
				/*if (adjacency_map.find(query_node2) == adjacency_map.end())
				{
					set<NODE_TYPE> temp_adj;
					temp_adj.insert((query_node1));
				}
				else
				{
					adjacency_map.find(query_node2)->second.insert(query_node1);
				}
				if (adjacency_map_reverse.find(query_node1) == adjacency_map_reverse.end())
				{
					set<NODE_TYPE> temp_adj;
					temp_adj.insert((query_node2));
				}
				else
				{
					adjacency_map_reverse.find(query_node1)->second.insert(query_node2);
				}*/



				// cout << "spd dfs lasts for " << double(end2 - start2) / CLOCKS_PER_SEC << endl;
				// cout << "hot points dfs lasts for " << double(end - start) / CLOCKS_PER_SEC - update_time << endl;
				// cout << "total order labels dfs lasts for " << double(end1 - end) / CLOCKS_PER_SEC << endl;
				// cout << "baseline dfs lasts for " << double(start - end2) / CLOCKS_PER_SEC << endl;




				result2.sort_by_string_order();
				result4.sort_by_string_order();
				result1.sort_by_string_order();
				if (output_result)
				{
					result1.write_to_file_append_edgeID("hp_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
					result2.write_to_file_append_edgeID("total_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
					result4.write_to_file_append_edgeID("dfs_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);
					result5.write_to_file_append_edgeID("bfs_pruned_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);
				}
				else
				{
					result1.write_to_file_append_edgeID_and_result_size("hp_result" + to_string(k) + "_" + dataset + "_" + to_string(threshold1), edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
					result2.write_to_file_append_edgeID_and_result_size("total_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);// +to_string(edgeIndex));
					result4.write_to_file_append_edgeID_and_result_size("dfs_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);
					result5.write_to_file_append_edgeID_and_result_size("bfs_pruned_result" + to_string(k) + "_" + dataset, edgeIndex, query_node1, query_node2);
				}
				pruned_lines += result5.path.size();
				baseline_lines += result1.path.size();
				out_performance << edgeIndex << " \t " << double(end2 - start2) / CLOCKS_PER_SEC << "\t" << double(end - start) / CLOCKS_PER_SEC - update_time << "\t " << double(end1 - end) / CLOCKS_PER_SEC << "\t" << double(start - end2) / CLOCKS_PER_SEC << "\t" << double(time_bf_pruned_path - end) / CLOCKS_PER_SEC << "\t" << double(time_bf_pruned_path_result5 - start2) / CLOCKS_PER_SEC << endl;

			}
			break;
	case 3:
//random querys
		cout << "generate random queries " << endl;
		out_performance << "edgeid \t " << algorithm << " query time " << endl;
		while (!inEdges.eof())
		{
			
			if (edgeIndex >= query_num)
			{
				break;
			}
			random_query_pair(query_node1, query_node2, node_num);
			//edgeIndex++;
			//// cout << edgeIndex << endl;
			//inEdges.getline(buffer, BUFFER_LENTH);
			//extractEdges(buffer, query_node2, query_node1);
			set<NODE_TYPE> c_path_set;
			vector<NODE_TYPE> c_path;
			//result3[0]
			paths result4, result5;
			paths al_result;
			if (edgeIndex == 0)
			{
				al_result.clear_file("result" + to_string(k) + "_" + algorithm + "_" + dataset);
			}
			start = clock();
			//choose algorithm
			int spd = bfs_return_spd(adjacency_list, adjacency_list_reverse, node_num, k, query_node1, query_node2);
			end = clock();
			if (spd == -1 || spd > k)
			{
				continue;
			}
			edgeIndex++;
			cout << edgeIndex << endl;
			//write the valid query edges to a specific file
			
			write_random_query_edges_withspd(dataset, k, query_node1, query_node2, edgeIndex,spd);
			


			//adjacency_list[query_node2].push_back(query_node1);//if the query is from node1 to node2, then the new edge is node2->node1 to from a circle
			//adjacency_list_reverse[query_node1].push_back(query_node2);

			
		}
		break;
	default:
		break;
	}

	



	

}  