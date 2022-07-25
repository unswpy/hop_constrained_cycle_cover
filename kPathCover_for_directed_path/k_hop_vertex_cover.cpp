#include "k_hop_vertex_cover.h"
set<NODE_TYPE> k_hop_vertex_cover(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num) {
	//copy adjacency_list
	vector<NODE_TYPE>* adjacency_list_temp = new vector<NODE_TYPE>[node_num+1];
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		vector<NODE_TYPE> temp(adjacency_list[i]);
		adjacency_list_temp[i].swap(temp);
	}
	set<NODE_TYPE> result;
	vector<hot_degree> hot_points;
	//double threhold = 0.05;//this means we find top 5% nodes in degree rank as hot-points(this part can change to k-core value, truess value or centrality)
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		NODE_TYPE temp_degree = adjacency_list[i].size()  + adjacency_list_reverse[i].size();// sum in and out degree
		if (temp_degree != 0)
		{
			hot_degree hd1(i, temp_degree);
			hot_points.push_back(hd1);
		}
	}
	stable_sort(hot_points.begin(), hot_points.end(), less<hot_degree>());
	for (NODE_TYPE j = 0; j < node_num; j++)
	{
		NODE_TYPE i = hot_points[j].node;
		//cout << i << endl;
		if (adjacency_list[i].size() == 0)
		{
			continue;
		}
		vector<NODE_TYPE> temp_k_path;
		vector<NODE_TYPE> c_path;
		if (find_k_path_from_node1(adjacency_list_temp, k, i, node_num, c_path, 0, temp_k_path))
		{//find a k_path, then push all the nodes in the k-path into result set and remove all the edges in and out from these nodes
			
			//// cout << "temp_k_path :";
			//for (auto iter1 = temp_k_path.begin(); iter1 != temp_k_path.end(); iter1++)
			//{
			//	// cout << *iter1 << " : ";
			//}
			//// cout << endl;
			for (auto cover_node = temp_k_path.begin(); cover_node != temp_k_path.end(); cover_node++)
			{
				result.insert(*cover_node);
				vector<NODE_TYPE> temp_v;
				adjacency_list_temp[*cover_node].swap(temp_v);//delete outgoing edges
				//start to delete in going edges
				for (auto in_node = adjacency_list_reverse[*cover_node].begin(); in_node != adjacency_list_reverse[*cover_node].end(); in_node++) {
					delete_adjacency_list_with_specific_edge(adjacency_list_temp, *cover_node, *in_node);
				}
			}
		}
	}
	delete[] adjacency_list_temp;
	return result;
}


bool node_neccessary(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result)//, set<NODE_TYPE>& block_nodes) 
{
	bool find_k_path_yet = false;	
	int cur_dis = 0;
	find_k_path_yet = find_k_path_from_node1(adjacency_list, k, node1, node_num, c_path, cur_dis, result);
	vector<NODE_TYPE> temp(result);
	c_path.swap(temp);
	cur_dis = c_path.size();
	if (find_k_path_yet) return true;
	else {
		for (auto iter2 = adjacency_list[node1].begin(); iter2 != adjacency_list[node1].end(); iter2++)
		{
			find_k_path_yet = find_k_path_from_node1(adjacency_list, k, *iter2, node_num, c_path, cur_dis, result);
			if (find_k_path_yet) return true;
		}
	}
	return false;
}


bool node_neccessary_enurate_all_paths(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result)//, set<NODE_TYPE>& block_nodes) 
{
	set<vector<NODE_TYPE>> temp_result;
	vector<NODE_TYPE> temp;
	temp.push_back(node1);
	bool find_k_path_yet = false;
	int cur_dis = 0;
	find_k_path_yet = find_k_path_from_node1_enum_paths(adjacency_list, k, node1, node_num, c_path, cur_dis, temp_result);
	temp_result.insert(temp);
	if (find_k_path_yet) {
		for (auto iter = temp_result.begin(); iter != temp_result.end(); iter++) {
			vector<NODE_TYPE> temp(*iter);
			if (temp.size()-1 >= k)
			{
				result.swap(temp);
				return true;
			}
		}
		return true;
	}
	else {
		for (auto iter_path = temp_result.begin(); iter_path != temp_result.end(); iter_path++)
		{
			vector<NODE_TYPE> temp(*iter_path);
			c_path.swap(temp);
			cur_dis = c_path.size();
			for (auto iter2 = adjacency_list_reverse[node1].begin(); iter2 != adjacency_list_reverse[node1].end(); iter2++)
			{
				find_k_path_yet = find_k_path_from_node1(adjacency_list_reverse, k, *iter2, node_num, c_path, cur_dis, result);
				if (find_k_path_yet)
				{
					vector<NODE_TYPE> temp2(*iter_path);
					vector<NODE_TYPE> temp3(result);
					for (int len = 0; len < temp3.size()- temp2.size(); len++)
					{
						result[len] = temp3[temp3.size() - len - 1];
					}
					for (int len = temp3.size() - temp2.size(); len < temp3.size(); len++)
					{
						result[len] = temp3[len];
					}
					return true;
				}
			}
		}
	}
	return false;
}


bool node_neccessary_enurate_all_paths_with_node_orders(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result,unordered_set<NODE_TYPE>& check_nodes, vector<NODE_TYPE>& node_orders)//, set<NODE_TYPE>& block_nodes) 
{
	set<vector<NODE_TYPE>> temp_result;
	bool find_k_path_yet = false;
	int cur_dis = 0;
	find_k_path_yet = find_k_path_from_node1_enum_paths_node_orders(adjacency_list, k, node1, node_num, c_path, cur_dis, temp_result, check_nodes,node_orders);

	if (find_k_path_yet) {
		for (auto iter = temp_result.begin(); iter != temp_result.end(); iter++) {
			vector<NODE_TYPE> temp(*iter);
			if (temp.size() - 1 >= k)
			{
				result.swap(temp);
				return true;
			}
		}
		return true;
	}
	else {
		for (auto iter_path = temp_result.begin(); iter_path != temp_result.end(); iter_path++)
		{
			vector<NODE_TYPE> temp(*iter_path);
			c_path.swap(temp);
			cur_dis = c_path.size();
			for (auto iter2 = adjacency_list[node1].begin(); iter2 != adjacency_list[node1].end(); iter2++)
			{
				if (check_nodes.find(*iter2) == check_nodes.end())
				{
					check_nodes.insert(*iter2);
					node_orders.push_back(*iter2);
				}
				find_k_path_yet = find_k_path_from_node1_node_orders(adjacency_list, k, *iter2, node_num, c_path, cur_dis, result,check_nodes,node_orders);
				if (find_k_path_yet) return true;
			}
		}
	}
	return false;
}


bool delete_adjacency_list_with_specific_edge(vector<NODE_TYPE>* adjacency_list, NODE_TYPE node1, NODE_TYPE node2) {
	bool result = false;
	for (auto iter = adjacency_list[node2].begin(); iter != adjacency_list[node2].end(); )
	{
		if (*iter == node1)
		{
			iter = adjacency_list[node2].erase(iter);
			result = true;
		}
		else
		{
			iter++;
		}
	}
	return result;
}


bool find_k_path_from_node1_enum_paths(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, set<vector<NODE_TYPE>>& result)//, set<NODE_TYPE>& block_nodes) 
{
	bool find_k_path_yet = false;
	//if node1 is already in c_path, then return false
	for (auto iter = c_path.begin(); iter != c_path.end(); iter++)
	{
		if (*iter == node1)
		{
			return false;
		}
	}
	c_path.push_back(node1);
	if (c_path.size() >= 2) {
		vector<NODE_TYPE> temp(c_path);
		result.insert(temp);
	}
	if (cur_distance >= k)
	{
		//vector<NODE_TYPE> temp(c_path);
		//result.push_back(temp);
		c_path.pop_back();
		return true;
	}
	for (auto iter = adjacency_list[node1].begin(); iter != adjacency_list[node1].end(); iter++)
	{
		find_k_path_yet = find_k_path_from_node1_enum_paths(adjacency_list, k, *iter, node_num, c_path, cur_distance + 1, result);
		if (find_k_path_yet)
		{
		break;
		}
	}
	c_path.pop_back();
	return find_k_path_yet;
}


bool find_k_path_from_node1_enum_paths_node_orders(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, set<vector<NODE_TYPE>>& result, unordered_set<NODE_TYPE>& check_nodes, vector<NODE_TYPE>& node_orders)//, set<NODE_TYPE>& block_nodes) 
{
	bool find_k_path_yet = false;
	//if node1 is already in c_path, then return false
	for (auto iter = c_path.begin(); iter != c_path.end(); iter++)
	{
		if (*iter == node1)
		{
			return false;
		}
	}
	c_path.push_back(node1);
	if (c_path.size() >= 2) {
		vector<NODE_TYPE> temp(c_path);
		result.insert(temp);
	}
	if (cur_distance >= k)
	{
		//vector<NODE_TYPE> temp(c_path);
		//result.push_back(temp);
		c_path.pop_back();
		return true;
	}
	for (auto iter = adjacency_list[node1].begin(); iter != adjacency_list[node1].end(); iter++)
	{
		if (check_nodes.find(*iter) == check_nodes.end())
		{
			check_nodes.insert(*iter);
			node_orders.push_back(*iter);
		}
		find_k_path_yet = find_k_path_from_node1_enum_paths_node_orders(adjacency_list, k, *iter, node_num, c_path, cur_distance + 1, result, check_nodes, node_orders);

		if (find_k_path_yet)
		{
			break;
		}
	}
	c_path.pop_back();
	return find_k_path_yet;
}




bool find_k_path_from_node1(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result)//, set<NODE_TYPE>& block_nodes) 
{
	bool find_k_path_yet = false;
	//if node1 is already in c_path, then return false
	for (auto iter = c_path.begin(); iter != c_path.end(); iter++)
	{
		if (*iter == node1)
		{
			vector<NODE_TYPE> temp(c_path);
			if (temp.size() > result.size())
			{
				result.swap(temp);
			}
			return false;
		}
	}
	c_path.push_back(node1);
	if (cur_distance >= k)
	{
		vector<NODE_TYPE> temp(c_path);
		result.swap(temp);
		c_path.pop_back();
		return true;
	}
	for (auto iter = adjacency_list[node1].begin(); iter != adjacency_list[node1].end(); iter++)
	{
		find_k_path_yet = find_k_path_from_node1(adjacency_list, k, *iter, node_num, c_path, cur_distance + 1, result);
		if (find_k_path_yet)
		{
			break;
		}
	}
	c_path.pop_back();
	return find_k_path_yet;
}


void find_bfs_order(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result, vector<bool>& check_nodes, vector<NODE_TYPE>& node_orders)//, set<NODE_TYPE>& block_nodes) 
{
	queue<NODE_TYPE> stk;
	stk.push(node);
	node_orders.push_back(node);
	check_nodes[node] = true;
	while (!stk.empty())
	{
		NODE_TYPE node1 = stk.front();
		stk.pop();
		
		
		//if node1 is already in c_path, then return false
		//auto iter = adjacency_list[node1].begin();
		for (auto iter = adjacency_list[node1].begin(); iter != adjacency_list[node1].end(); iter++)
		{

			if (!check_nodes[*iter])
			{

				node_orders.push_back(node1);

				stk.push(*iter);
				check_nodes[*iter] = true;
				//break;
			}
		}
		//if (iter == adjacency_list[node1].end())
		//{
		//	stk.pop();
		//	node_orders.push_back(node1);
		//}
	}
	return;
}


void find_dfs_order(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result, vector<bool>& check_nodes, vector<NODE_TYPE>& node_orders)//, set<NODE_TYPE>& block_nodes) 
{
	stack<NODE_TYPE> stk;
	stk.push(node);
	check_nodes[node] = true;
	while (!stk.empty())
	{
		NODE_TYPE node1 = stk.top();
		//stk.pop();
	
			check_nodes[node1] = true;
			//if node1 is already in c_path, then return false
			auto iter = adjacency_list[node1].begin();
			for (iter = adjacency_list[node1].begin(); iter != adjacency_list[node1].end(); iter++)
			{

				if (!check_nodes[*iter])
				{
					stk.push(*iter);
					break;
				}
			}
			if (iter == adjacency_list[node1].end())
			{
				stk.pop();
				node_orders.push_back(node1);
			}
	}
	return ;
}


void find_dfs_order_pre(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result, vector<bool>& check_nodes, vector<NODE_TYPE>& node_orders)//, set<NODE_TYPE>& block_nodes) 
{
	stack<NODE_TYPE> stk;
	stk.push(node);
	//check_nodes[node] = true;
	while (!stk.empty())
	{
		NODE_TYPE node1 = stk.top();
		//stk.pop();
		if (!check_nodes[node1])
		{
			node_orders.push_back(node1);
		}
		check_nodes[node1] = true;
		

		//if node1 is already in c_path, then return false
		auto iter = adjacency_list[node1].begin();
		for (iter = adjacency_list[node1].begin(); iter != adjacency_list[node1].end(); iter++)
		{

			if (!check_nodes[*iter])
			{
				stk.push(*iter);
				//node_orders.push_back(node1);
				break;
			}
		}
		if (iter == adjacency_list[node1].end())
		{
			stk.pop();
		}
	}
	return;
}



bool find_k_path_from_node1_node_orders(vector<NODE_TYPE>* adjacency_list, int k, NODE_TYPE node1, NODE_TYPE node_num, vector<NODE_TYPE>& c_path, int cur_distance, vector<NODE_TYPE>& result, unordered_set<NODE_TYPE>& check_nodes, vector<NODE_TYPE>& node_orders)//, set<NODE_TYPE>& block_nodes) 
{
	bool find_k_path_yet = false;
	//if node1 is already in c_path, then return false
	for (auto iter = c_path.begin(); iter != c_path.end(); iter++)
	{
		if (*iter == node1)
		{
			vector<NODE_TYPE> temp(c_path);
			if (temp.size() > result.size())
			{
				result.swap(temp);
			}
			return false;
		}
	}
	c_path.push_back(node1);
	if (cur_distance >= k)
	{
		vector<NODE_TYPE> temp(c_path);
		result.swap(temp);
		c_path.pop_back();
		return true;
	}
	for (auto iter = adjacency_list[node1].begin(); iter != adjacency_list[node1].end(); iter++)
	{
		if (check_nodes.find(*iter) == check_nodes.end())
		{
			check_nodes.insert(*iter);
			node_orders.push_back(*iter);
		}
		find_k_path_yet = find_k_path_from_node1_node_orders(adjacency_list, k, *iter, node_num, c_path, cur_distance + 1, result,check_nodes,node_orders);

		if (find_k_path_yet)
		{
			break;
		}
	}
	c_path.pop_back();
	return find_k_path_yet;
}



NODE_TYPE find_result_node_in_single_k_path(vector<NODE_TYPE> temp_k_path, int* find_times)
{
	NODE_TYPE temp_index = temp_k_path.size() / 2;
	NODE_TYPE temp_max = find_times[temp_k_path[temp_index]];
	NODE_TYPE temp_max_index = temp_index;
	NODE_TYPE left, right;
	left = temp_index - 1;
	right = temp_index;
	while (left >= 0 )
	{
		if (left >= 0)
		{
			if (find_times[temp_k_path[left]] > temp_max) {
				temp_max_index = left;
				temp_max = find_times[temp_k_path[left]];
			}
			left--;
		}
	}
	while (right < temp_k_path.size())
	{
		if (right < temp_k_path.size())
		{
			if (find_times[temp_k_path[right]] > temp_max) {
				temp_max_index = right;
				temp_max = find_times[temp_k_path[right]];
			}
			right++;
		}
	}
	return temp_k_path[temp_max_index];
}

unordered_set<NODE_TYPE> k_hop_vertex_cover_heuristic(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num)
{
	node_num++;
	//copy adjacency_list
	int* find_times = new int[node_num+1];
	//cout << node_num << endl;
	vector<NODE_TYPE>* adjacency_list_temp = new vector<NODE_TYPE>[node_num + 1];
	vector<NODE_TYPE>* adjacency_list_temp_reverse = new vector<NODE_TYPE>[node_num + 1];
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		find_times[i] = 0;
		vector<NODE_TYPE> temp(adjacency_list[i]);
		adjacency_list_temp[i].swap(temp);
		vector<NODE_TYPE> temp2(adjacency_list_reverse[i]);
		adjacency_list_temp_reverse[i].swap(temp2);
	}
	unordered_set<NODE_TYPE> result;
	vector<hot_degree> hot_points;
	//double threhold = 0.05;//this means we find top 5% nodes in degree rank as hot-points(this part can change to k-core value, truess value or centrality)
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		NODE_TYPE temp_degree = adjacency_list[i].size()  + adjacency_list_reverse[i].size();// sum in and out degree
		//if (temp_degree != 0)
		{
			hot_degree hd1(i, temp_degree);
			hot_points.push_back(hd1);
		}
	}
	stable_sort(hot_points.begin(), hot_points.end(), less<hot_degree>());
	unordered_set<NODE_TYPE> second_iter;

	vector<bool> check_nodes;
	check_nodes.resize(node_num + 1);
	for (int i = 0; i <= node_num; i++)
	{
		check_nodes[i] = false;
	}
	vector<NODE_TYPE> temp_k_path1;
	vector<NODE_TYPE> c_path1;
	vector<NODE_TYPE> node_orders;
	//sort all the vertex by the dfs order
	for (int i = 0; i < node_num; i++)
	{
		if (check_nodes[i])
			continue;
		find_dfs_order_pre(adjacency_list, k, i, node_num, c_path1, 0, temp_k_path1, check_nodes, node_orders);
	}

	vector<NODE_TYPE> temp_node_orders;

	cout << node_orders.size() << " order\t";

	cout << endl;
	//for (auto node_ex = node_orders.begin(); node_ex != node_orders.end(); node_ex++)
	//{


	for (int j = 0; j < node_num; j++)
		//for (auto node_ex = node_orders.rbegin(); node_ex != node_orders.rend(); node_ex++)
	{
		NODE_TYPE i = hot_points[j].node; //*node_ex;//hot_points[i].node;
		//NODE_TYPE i = *node_ex;
		//NODE_TYPE i =  hot_points[j].node; //*node_ex;//= hot_points[j].node;
		if (result.find(i) != result.end()) continue;
		if (i == 328755)
		{
			cout << 328755 << endl;
		}
		//if (hot_points[j].degree == 0) continue;
		//if (find_times[i] != 0) 
		//{
		//	second_iter.insert(i);
		//	continue;
		//}
		//if(i % 300000 == 0) cout << i << endl;
		//if (adjacency_list[i].size() == 0)
		//{
		//	continue;
		//}
		vector<NODE_TYPE> temp_k_path;
		vector<NODE_TYPE> c_path;
		if (node_neccessary_enurate_all_paths(adjacency_list_temp, adjacency_list_temp_reverse, k, i, node_num, c_path, 0, temp_k_path))//find_k_path_from_node1(adjacency_list_temp, k, i, node_num, c_path, 0, temp_k_path))
		{//find a k_path, then push all the nodes in the k-path into result set and remove all the edges in and out from these nodes

		 //// cout << "temp_k_path :";
		 //for (auto iter1 = temp_k_path.begin(); iter1 != temp_k_path.end(); iter1++)
		 //{
		 //	// cout << *iter1 << " : ";
		 //}
		 //// cout << endl;

			for (auto cover_node = temp_k_path.begin(); cover_node != temp_k_path.end(); cover_node++)
			{
				find_times[*cover_node]++;
			}
			NODE_TYPE result_node = 0;

			result_node = find_result_node_in_single_k_path(temp_k_path, find_times);
			
			result.insert(result_node);
			vector<NODE_TYPE> temp_v;
			adjacency_list_temp[result_node].swap(temp_v);//delete outgoing edges
														  //start to delete in going edges
			vector<NODE_TYPE> temp_v2;
			adjacency_list_temp_reverse[result_node].swap(temp_v2);
			for (auto in_node = adjacency_list_reverse[result_node].begin(); in_node != adjacency_list_reverse[result_node].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp, result_node, *in_node);
			}
			for (auto in_node = adjacency_list[result_node].begin(); in_node != adjacency_list[result_node].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp_reverse, result_node, *in_node);
			}


			//node_ex--;
			j--;
		}
	}
	//second iter
	/*
	for (auto iter = second_iter.begin(); iter != second_iter.end(); )
	{
		NODE_TYPE i = *iter;
		if (result.find(i) != result.end()) continue;
		if (i == 328755)
		{
			cout << 328755 << endl;
			vector<int> w = { 1312646,328760,1319119,332238 };
			for (auto num = w.begin(); num != w.end(); num++)
			{
				if (result.find(*num) != result.end())
				{
					cout << *num << "is in the result set" << endl;
				}
			}
			cout << "check 328755 done" << endl;
		}
		if (i % 300000 == 0) cout << i << endl;
		//if (adjacency_list[i].size() == 0)
		//{
		//	continue;
		//}
		vector<NODE_TYPE> temp_k_path;
		vector<NODE_TYPE> c_path;
		//if (find_k_path_from_node1(adjacency_list_temp, k, i, node_num, c_path, 0, temp_k_path))
		if (node_neccessary_enurate_all_paths(adjacency_list_temp, k, i, node_num, c_path, 0, temp_k_path))
		{//find a k_path, then push all the nodes in the k-path into result set and remove all the edges in and out from these nodes

		 //// cout << "temp_k_path :";
		 //for (auto iter1 = temp_k_path.begin(); iter1 != temp_k_path.end(); iter1++)
		 //{
		 //	// cout << *iter1 << " : ";
		 //}
		 //// cout << endl;

			for (auto cover_node = temp_k_path.begin(); cover_node != temp_k_path.end(); cover_node++)
			{
				find_times[*cover_node]++;
			}
			NODE_TYPE reuslt_node = find_result_node_in_single_k_path(temp_k_path, find_times);
			result.insert(reuslt_node);
			vector<NODE_TYPE> temp_v;
			adjacency_list_temp[reuslt_node].swap(temp_v);//delete outgoing edges
														  //start to delete in going edges
			for (auto in_node = adjacency_list_reverse[reuslt_node].begin(); in_node != adjacency_list_reverse[reuslt_node].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp, reuslt_node, *in_node);
			}

			
		}
		else {
			iter++;
		}
	}
	*/
	// try to get the minimal one
	/*
	for (int w = 0; w < 1; w++)
	{
		for (auto iter = result.begin(); iter != result.end(); )
		{
			if (*iter == 1312646)
			{
				cout << 1312646 << endl;
			}
			//try to delete *iter
			for (auto out_node = adjacency_list_reverse[*iter].begin(); out_node != adjacency_list_reverse[*iter].end(); out_node++) {
				if (result.find(*out_node) != result.end()) continue;
				//ector<NODE_TYPE> temp_nv(adjacency_list[*iter]);
				adjacency_list_temp[*iter].push_back(*out_node);
			}
			for (auto in_node = adjacency_list_reverse[*iter].begin(); in_node != adjacency_list_reverse[*iter].end(); in_node++) {
				if (result.find(*in_node) != result.end()) continue;
				adjacency_list_temp[*in_node].push_back(*iter);
			}




			vector<NODE_TYPE> temp_k_path;
			vector<NODE_TYPE> c_path;
			if (!find_k_path_from_node1(adjacency_list_temp, k, *iter, node_num, c_path, 0, temp_k_path))
			{// need to explore two directions
				iter = result.erase(iter);
			}
			else {
				vector<NODE_TYPE> temp_nv2;
				adjacency_list_temp[*iter].swap(temp_nv2);
				for (auto in_node = adjacency_list_reverse[*iter].begin(); in_node != adjacency_list_reverse[*iter].end(); in_node++) {
					delete_adjacency_list_with_specific_edge(adjacency_list_temp, *iter, *in_node);
				}
				iter++;
			}
		}


		for (NODE_TYPE j = 0; j < node_num; j++)
		{
			NODE_TYPE i = hot_points[j].node;
			if (i == 1312646)
			{
				cout << 1312646 << endl;
			}
			if (result.find(i) != result.end())
			{
				continue;
			}
			//if (find_times[i] != 0)
			//{
			//	second_iter.insert(i);
			//	continue;
			//}
			if (i % 300000 == 0) cout << i << endl;
			if (adjacency_list[i].size() == 0)
			{
				continue;
			}
			vector<NODE_TYPE> temp_k_path;
			vector<NODE_TYPE> c_path;
			if (find_k_path_from_node1(adjacency_list_temp, k, i, node_num, c_path, 0, temp_k_path))
			{//find a k_path, then push all the nodes in the k-path into result set and remove all the edges in and out from these nodes

			 //// cout << "temp_k_path :";
			 //for (auto iter1 = temp_k_path.begin(); iter1 != temp_k_path.end(); iter1++)
			 //{
			 //	// cout << *iter1 << " : ";
			 //}
			 //// cout << endl;

				for (auto cover_node = temp_k_path.begin(); cover_node != temp_k_path.end(); cover_node++)
				{
					find_times[*cover_node]++;
				}
				NODE_TYPE result_node = 0;
				if (j <= 0)
				{
					result_node = i;
				}
				else {
					result_node =  find_result_node_in_single_k_path(temp_k_path, find_times);
				}
				result.insert(result_node);
				vector<NODE_TYPE> temp_v;
				adjacency_list_temp[result_node].swap(temp_v);//delete outgoing edges
															  //start to delete in going edges
				for (auto in_node = adjacency_list_reverse[result_node].begin(); in_node != adjacency_list_reverse[result_node].end(); in_node++) {
					delete_adjacency_list_with_specific_edge(adjacency_list_temp, result_node, *in_node);
				}


			}
		}
	}
	*/
	//minimal
	//for (auto iter = node_orders.begin(); iter != node_orders.end(); iter++)
	//{
	
	//for(int i = 0; i < node_num;i++)
	for (auto node_ex = node_orders.rbegin(); node_ex != node_orders.rend(); node_ex++)
	{
		NODE_TYPE j = *node_ex; //hot_points[i].node; //*node_ex;//hot_points[i].node;
		if (result.find(j) == result.end()) continue;
		//try to delete *iter
		for (auto out_node = adjacency_list[j].begin(); out_node != adjacency_list[j].end(); out_node++) {
			if (result.find(*out_node) != result.end()) continue;
			//ector<NODE_TYPE> temp_nv(adjacency_list[*iter]);
			adjacency_list_temp[j].push_back(*out_node);
			adjacency_list_temp_reverse[*out_node].push_back(j);
		}
		for (auto in_node = adjacency_list_reverse[j].begin(); in_node != adjacency_list_reverse[j].end(); in_node++) {
			if (result.find(*in_node) != result.end()) continue;
			adjacency_list_temp_reverse[j].push_back(*in_node);
			adjacency_list_temp[*in_node].push_back(j);
		}


		vector<NODE_TYPE> temp_k_path;
		vector<NODE_TYPE> c_path;
		if (!node_neccessary_enurate_all_paths(adjacency_list_temp, adjacency_list_temp_reverse, k, j, node_num, c_path, 0, temp_k_path))
		{// need to explore two directions
			result.erase(result.find(j));
		}
		else {
			if (j == 263555)
			{
				for (auto temp = temp_k_path.begin(); temp != temp_k_path.end(); temp++)
				{
					cout << *temp << "\t";
				}
				cout << endl;
			}
			vector<NODE_TYPE> temp_nv2;
			adjacency_list_temp[j].swap(temp_nv2);
			vector<NODE_TYPE> temp_nv3;
			adjacency_list_temp_reverse[j].swap(temp_nv3);
			for (auto in_node = adjacency_list_reverse[j].begin(); in_node != adjacency_list_reverse[j].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp, j, *in_node);
			}
			for (auto in_node = adjacency_list[j].begin(); in_node != adjacency_list[j].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp_reverse, j, *in_node);
			}

			//iter++;
		}
	}
	
	cout << "iter1\t" << result.size() << endl;
	delete[] adjacency_list_temp_reverse;
	delete[] adjacency_list_temp;
	delete[] find_times;
	return result;
}

bool judge_cover_correctness(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num, unordered_set<NODE_TYPE>& test_cover)
{
	unordered_set<NODE_TYPE> result(test_cover);
	vector<NODE_TYPE>* adjacency_list_temp = new vector<NODE_TYPE>[node_num + 1];
	vector<NODE_TYPE>* adjacency_list_temp_reverse = new vector<NODE_TYPE>[node_num + 1];
	for (int i = 0; i < node_num; i++)
	{
		if (i == 1)
		{
			cout << 1 << "endl";
		}
		if (result.find(i) != result.end()) continue;
		for (auto out_node = adjacency_list[i].begin(); out_node != adjacency_list[i].end(); out_node++) {
			if (result.find(*out_node) != result.end()) continue;
			//ector<NODE_TYPE> temp_nv(adjacency_list[*iter]);
			adjacency_list_temp[i].push_back(*out_node);
			adjacency_list_temp_reverse[*out_node].push_back(i);
		}
		for (auto in_node = adjacency_list_reverse[i].begin(); in_node != adjacency_list_reverse[i].end(); in_node++) {
			if (result.find(*in_node) != result.end()) continue;
			adjacency_list_temp_reverse[i].push_back(*in_node);
			adjacency_list_temp[*in_node].push_back(i);

		}

	}
	
	cout << "before minial test iter1\t" << result.size() << endl;
	unordered_set<NODE_TYPE> result2(result);
	//minimal
	
	//for (auto iter = result.begin(); iter != result.end(); )
	//{
	unordered_set<NODE_TYPE> nums;
	nums.insert(0);
	nums.insert(3143930);
	nums.insert(3143888);
	nums.insert(3143887);

	for (int j = 0; j < node_num; j++)
	{

		//NODE_TYPE j = hot_points[i].node;
		if (result.find(j) == result.end()) continue;
		if (nums.find(j) != nums.end())
		{
			cout << j << endl;
		}
		//try to delete *iter
		for (auto out_node = adjacency_list[j].begin(); out_node != adjacency_list[j].end(); out_node++) {
			if (result.find(*out_node) != result.end()) continue;
			//ector<NODE_TYPE> temp_nv(adjacency_list[*iter]);
			adjacency_list_temp[j].push_back(*out_node);
			adjacency_list_reverse[*out_node].push_back(j);
		}
		for (auto in_node = adjacency_list_reverse[j].begin(); in_node != adjacency_list_reverse[j].end(); in_node++) {
			if (result.find(*in_node) != result.end()) continue;
			adjacency_list_temp_reverse[j].push_back(*in_node);
			adjacency_list[*in_node].push_back(j);
		}


		vector<NODE_TYPE> temp_k_path;
		vector<NODE_TYPE> c_path;
		if (!node_neccessary_enurate_all_paths(adjacency_list_temp, adjacency_list_temp_reverse, k, j, node_num, c_path, 0, temp_k_path))
		{// need to explore two directions
			result.erase(result.find(j));
		}
		else {
			if (j == 263555)
			{
				for (auto temp = temp_k_path.begin(); temp != temp_k_path.end(); temp++)
				{
					cout << *temp << "\t";
				}
				cout << endl;
			}
			vector<NODE_TYPE> temp_nv2;
			adjacency_list_temp[j].swap(temp_nv2);
			vector<NODE_TYPE> temp_nv3;
			adjacency_list_temp_reverse[j].swap(temp_nv3);
			for (auto in_node = adjacency_list_reverse[j].begin(); in_node != adjacency_list_reverse[j].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp, j, *in_node);
			}
			for (auto in_node = adjacency_list[j].begin(); in_node != adjacency_list[j].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp_reverse, j, *in_node);
			}

			//iter++;
		}
	}
	
	if (result2.size() != result.size())
	{
		cout << "the cover is not the minimal" << endl;
		cout << "size of minial test" << result.size() << endl;
		//return false;
	}
	
	for (int i = 0; i < node_num; i++ )
	{
		if (result.find(i) != result.end()) continue;
		//try to delete *iter
		vector<NODE_TYPE> temp_k_path;
		vector<NODE_TYPE> c_path;
		if (i == 1 || i == 500600 || i == 143280 || i == 320348)
		{
			cout << 1 << endl;
			vector<int> w = { 1,500600,143280,320348 };
			for (auto num = w.begin(); num != w.end(); num++)
			{
				if (result.find(*num) != result.end())
				{
					cout << *num << "is in the result set" << endl;
				}
			}
			cout << "check 328755 done" << endl;
			//if(result.find(32))
		}
		if (node_neccessary_enurate_all_paths(adjacency_list_temp, adjacency_list_temp_reverse, k, i, node_num, c_path, 0, temp_k_path))
		{// need to explore two directions
		 //iter = result.erase(iter);
			cout << i << endl;
			for (auto iter = temp_k_path.begin(); iter != temp_k_path.end(); iter++)
			{
				cout << *iter << "\t";
			}
			cout << endl;
			cout << "node is " << i << endl;
			cout << "iter2\t" << result.size() << endl;
			cout << "not a valid k path cover" << endl;
			return false;
		}
	}
	cout << "iter2\t" << result.size() << endl;

	delete[] adjacency_list_temp_reverse;
	delete[] adjacency_list_temp;
	return true;
}




unordered_set<NODE_TYPE> k_hop_vertex_cover_node_necessary_dfs_order(vector<NODE_TYPE>* adjacency_list, vector<NODE_TYPE>* adjacency_list_reverse, int k, NODE_TYPE node_num)
{
	//copy adjacency_list
	//int* find_times = new int[node_num + 1];
	vector<NODE_TYPE>* adjacency_list_temp = new vector<NODE_TYPE>[node_num + 1];
	vector<NODE_TYPE>* adjacency_list_temp_reverse = new vector<NODE_TYPE>[node_num + 1];
	//for (NODE_TYPE i = 0; i < node_num; i++)
	//{
	//vector<NODE_TYPE> temp(adjacency_list[i]);
	//adjacency_list_temp[i].swap(temp);
	//}
	unordered_set<NODE_TYPE> result;
	vector<hot_degree> hot_points;
	//double threhold = 0.05;//this means we find top 5% nodes in degree rank as hot-points(this part can change to k-core value, truess value or centrality)
	for (NODE_TYPE i = 0; i < node_num; i++)
	{
		NODE_TYPE temp_degree = adjacency_list[i].size() + adjacency_list_reverse[i].size();// sum in and out degree
																							//if (temp_degree != 0)
		{
			hot_degree hd1(i, temp_degree);
			hot_points.push_back(hd1);
		}
	}
	stable_sort(hot_points.begin(), hot_points.end(), less<hot_degree>());
	set<NODE_TYPE> second_iter;
	//unordered_set<NODE_TYPE> node_orders;
	vector<bool> tested;
	tested.resize(node_num);
	for (int i = 0; i < node_num; i++)
	{
		tested[i] = false;
	}
	for (NODE_TYPE j = 0; j < node_num; j++)
	{

		NODE_TYPE i = hot_points[j].node;
		//if (i == 3) continue;
		result.insert(i);
	}
	//result.insert(3);
	//second iter
	vector<NODE_TYPE> node_orders;
	vector<bool> check_nodes;
	check_nodes.resize(node_num + 1);
	for (int i = 0; i < node_num + 1; i++)
	{
		check_nodes[i] = false;
	}
	unordered_set<NODE_TYPE> result_temp(result);

	int sample_vertex = 0;
	vector<NODE_TYPE> temp_k_path1;
	vector<NODE_TYPE> c_path1;
	//sort all the vertex by the dfs order
	for (int j = 0; j < node_num; j++)
	{
		NODE_TYPE i = hot_points[j].node;
		if (check_nodes[i])
			continue;
		find_dfs_order(adjacency_list, k, i, node_num, c_path1, 0, temp_k_path1, check_nodes, node_orders);
	}

	vector<NODE_TYPE> temp_node_orders;

	cout << node_orders.size() << " order\t";

	cout << endl;
	for (auto node_ex = node_orders.begin(); node_ex != node_orders.end(); node_ex++)
	{
		if (*node_ex == 98 || *node_ex == 466991 || *node_ex == 491346 || *node_ex == 322885 )// || *node_ex == 500600 || *node_ex == 143280 || *node_ex == 320346 )
		{
			cout << *node_ex << endl;
		}
		if (tested[*node_ex]) continue;
		tested[*node_ex] = true;
		for (auto out_node = adjacency_list[*node_ex].begin(); out_node != adjacency_list[*node_ex].end(); out_node++) {
			if (result.find(*out_node) != result.end()) continue;
			//ector<NODE_TYPE> temp_nv(adjacency_list[*iter]);
			adjacency_list_temp[*node_ex].push_back(*out_node);
			adjacency_list_temp_reverse[*out_node].push_back(*node_ex);
		}
		for (auto in_node = adjacency_list_reverse[*node_ex].begin(); in_node != adjacency_list_reverse[*node_ex].end(); in_node++) {
			if (result.find(*in_node) != result.end()) continue;
			adjacency_list_temp_reverse[*node_ex].push_back(*in_node);
			adjacency_list_temp[*in_node].push_back(*node_ex);
		}


		vector<NODE_TYPE> temp_k_path;
		vector<NODE_TYPE> c_path;
		if (!node_neccessary_enurate_all_paths(adjacency_list_temp, adjacency_list_temp_reverse, k, *node_ex, node_num, c_path, 0, temp_k_path))
		{// need to explore two directions
		 //iter = result.erase(iter);
			result.erase(result.find(*node_ex));
		}
		else {
			vector<NODE_TYPE> temp_nv2;
			adjacency_list_temp[*node_ex].swap(temp_nv2);
			vector<NODE_TYPE> temp_nv3;
			adjacency_list_temp_reverse[*node_ex].swap(temp_nv3);
			for (auto in_node = adjacency_list_reverse[*node_ex].begin(); in_node != adjacency_list_reverse[*node_ex].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp, *node_ex, *in_node);
			}
			for (auto in_node = adjacency_list[*node_ex].begin(); in_node != adjacency_list[*node_ex].end(); in_node++) {
				delete_adjacency_list_with_specific_edge(adjacency_list_temp_reverse, *node_ex, *in_node);
			}
		}
		



	}


	delete[] adjacency_list_temp;
	return result;
}



void test_k_hop_cover(int k, char* dataset, char* algorithm)
{//algorithms includes 
	k = k -1;
	//int option = 0;// 0 means use dynamic files to test
	//double threshold1 = 0.0001;//hot points threshold
	//double threshold2 = 0.0010;//label nodes threshold
	//int k = 5;
	// 1 means test with given query_node1 and query_node2
	ofstream test;
	//test.open("Amazon_result\\hp_result" + to_string(testi));
	//test << "test " << endl;


	ifstream count_file;
	ofstream output;
	ifstream inEdges;
	string outfile = dataset;
	outfile += "_result"+ to_string(k);
	NODE_TYPE  node_num;// = 548552; //121790;// 548552; //1007518272;//548552;//1234944765;
						//NODE_TYPE node_num = 91306;
	string count_name = dataset;
	count_file.open(count_name);
	output.open(outfile);
	char tmp[BUFFER_LENTH];

	if (!count_file.is_open())
	{
		cout << "Error opening file : " << dataset << endl;
	}
	//configure the values of n and m 
	node_num = 0;
	int i = 0;
	NODE_TYPE x, y;
	bool roadUSA = false;
	bool pound = false;
	while (!count_file.eof())
	{
		NODE_TYPE x, y;
		i++;
		count_file.getline(tmp, BUFFER_LENTH);
		if (i == 1 && tmp[0] == 'c')
		{
			roadUSA = true;
		}
		else if (i == 1 && (tmp[0] == '#' || (tmp[0] >= '0' && tmp[0] <= '9')))
		{
			pound = true;
		}
		if (roadUSA)
		{
			if (tmp[0] != 'a') continue;
			extractEdgesSkipNodes(tmp, x, y, 1);
		}
		else if (pound)
		{
			if (tmp[0] == '#') continue;
			extractEdges(tmp, x, y);
		}
		else {

		}
		if (i <= 10)
		{
			cout << x << " : " << y << endl;
		}
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

	

	vector<int> in_degrees, out_degrees;// we need to get in and out degrees
	vector<NODE_TYPE >* adjacency_list = new vector<NODE_TYPE >[node_num + 1];
	vector<NODE_TYPE >* adjacency_list_reverse = new vector<NODE_TYPE >[node_num + 1];
	vector<NODE_TYPE >* adjacency_list_double = new vector<NODE_TYPE>[node_num + 1];//record the double directed adjacency for the undirected graph load 

	string temp_static(dataset);

	load_graph_data(temp_static, adjacency_list, adjacency_list_reverse);
	NODE_TYPE node1,cur_distance;
	vector<NODE_TYPE> c_path,result;
	//node_neccessary(adjacency_list, k, node1, node_num, c_path, cur_distance, result);
	unordered_set<NODE_TYPE> result1, result3;
	unordered_set<NODE_TYPE> result2;
	clock_t t1,t2,t3;
	t1 = clock();
	result1 = k_hop_vertex_cover_heuristic(adjacency_list, adjacency_list_reverse, k, node_num);
	//result1 = k_hop_vertex_cover_node_necessary_dfs_pre_order(adjacency_list, adjacency_list_reverse, k, node_num);
	//result1 = k_hop_vertex_cover_node_necessary_bfs_order(adjacency_list, adjacency_list_reverse, k, node_num);
	t2 = clock();
	cout << "baseline:\t" << result2.size() << "\t ours : \t" << result1.size() << endl;
	cout << "runtime:ours\t" << double(t2 - t1) / CLOCKS_PER_SEC << endl;
	cout << "ours done" << endl;
	output << "baseline:\t" << result2.size() << "\t ours : \t" << result1.size() << endl;
	output << "runtime:ours\t" << double(t2 - t1) / CLOCKS_PER_SEC << endl;
	output << "ours done" << endl;
	result2 = k_hop_vertex_cover_node_necessary_dfs_order(adjacency_list, adjacency_list_reverse, k, node_num);
	t3 = clock();
	cout << "baseline:\t" << result2.size() <<"\t ours : \t" << result1.size() << endl;
	cout << "runtime:ours\t" << double(t2 - t1) / CLOCKS_PER_SEC << "\t baseline : \t" << double(t3 - t2) / CLOCKS_PER_SEC << endl;
	output << "baseline:\t" << result2.size() << "\t ours : \t" << result1.size() << endl;
	output << "runtime:ours\t" << double(t2 - t1) / CLOCKS_PER_SEC << "\t basline : \t" << double(t3 - t2) / CLOCKS_PER_SEC << endl;

	bool test1 =  judge_cover_correctness(adjacency_list, adjacency_list_reverse, k, node_num, result1);
	vector<int> w = { 7,1048583,334,1048909,1048885,1048884 };
	for (auto num = w.begin(); num != w.end(); num++)
	{
		if (result1.find(*num) != result1.end())
		{
			cout << *num << "is in the result1 set" << endl;
		}
	}
	cout << "result1 check 2105 done" << endl;



	vector<int> w2 = { 98,466991,491346,322885 };
	for (auto num = w2.begin(); num != w2.end(); num++)
	{
		if (result2.find(*num) != result2.end())
		{
			cout << *num << "is in the result2 set 233" << endl;
		}
	}
	cout << "result2 check 1 done" << endl;
	bool test2 = judge_cover_correctness(adjacency_list, adjacency_list_reverse, k, node_num, result2);
	if (!test1)
	{
		cout << "ours are wrong" << dataset<< endl;
	}
	else {
		cout << "ours are correct" << endl;
	}
	if (!test2)
	{
		cout << "baseline is wrong" << endl;
	}
	else {
		cout << "baseline is correct" << endl;
	}
	if (result1.size() <= 10 && result2.size() <= 10)
	{
		for (auto iter = result1.begin(); iter != result1.end(); iter++)
		{
			cout << *iter << "\t";
		}
		cout << "end of result1 " << endl;

		for (auto iter = result2.begin(); iter != result2.end(); iter++)
		{
			cout << *iter << "\t";
		}
		cout << "end of result2 " << endl;
	}
	system("pause");
	output.close();
}
