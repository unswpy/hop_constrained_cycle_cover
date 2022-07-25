#include "k_hop_vertex_cover.h"
#include "k_reach_total_order.h"
int main()
{
	////join_drop_longpaths_and_repeat_nodes_and_short_paths_stop_ndoeswithdis
	//vector<long long> a;
	//a.push_back(1);
	//a.push_back(2);
	//a.push_back(3);
	//a.push_back(4);
	//vector<long long> b;
	//b.push_back(4);
	//b.push_back(5);
	//b.push_back(6);
	//paths pa;
	//paths pb;
	//pa.push_back(a);
	//pb.push_back(b);
	//set<NODE_TYPE> stop_nodes;
	////stop_nodes.insert(5);
	//paths result;
	//pa.join_drop_longpaths_and_repeat_nodes_and_short_paths_stop_ndoeswithdis(pb, 15, result, 4, stop_nodes, 1, 1);
	//result.output();
	//system("pause");


	//vector <NODE_TYPE > v;
	//for(auto iter_v = a.begin()+1; iter_v != a.end()-1; iter_v++)
	//{
	//	// cout << *iter_v << endl;
	//}
	//system("pause");
	//bug with dfs_with_block_induced_spd_ini
	char datasets[255][100] = {"advogato.edge", "Amazon.txt","twitter_social","gplus_combined.txtreoder", "Slashdot0902.txtreoder",
		"soc-Epinions1.txtreoder", "soc-LiveJournal1.txtreoder","soc-pokec-relationships.txtreoder","WikiTalk.txtreoder" };
	bool output_result = false;
	//0 means 
	char dataset[255] = "cit-Patents.txtreoder";
	int k = 6;

	// 1 means test with given query_node1 and query_node2
	//test_k_hop_cover(5, "Amazon.txt.static", "ours");
	test_k_hop_cover(16, "CAL", "ours");
	//test_k_hop_cover(4, "test", "ours");
	//system("pause");
	return 0;

}

