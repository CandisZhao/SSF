#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <hash_set>
#include <hash_map>
#include <algorithm>
#include <iostream>
using namespace std;


class GraphCluster
{
public:
	
	long double alpha = 0.01;    // 参数
	string in_file;
	string out_file;

	vector<long double> f;


	int n;                        // 图节点数目
	int D;
	int DS;
	hash_map<string, int> string2int;
	hash_map<int, string> int2string;

	vector< vector <int> > graph;   // 邻接表
	vector<int> sub_graph_node;
	vector<hash_set<int>> clusters;     // 子图节点

	vector<int> ki;                 // 节点度
	vector<int> ki_in;              // 节点与子图相邻度
	vector<int> ki_out;             // 节点度-节点与子图相邻度

	vector< vector <int> > S;
	vector< vector <int> > F;
	vector<bool> visited;           // 已遍历的节点


	GraphCluster(long double alp, string input, string output);
	void initialize();              // 初始化
	void initialize2();

	void SSF();   // 主函数
	void SSF_old();   // 主函数
	vector<int> search_one_complex();  
	int max_degree_unvisited();    
	void Setcover();                // 后处理

	bool is_same_graph(hash_set<int> &graph1, hash_set<int>&graph2);
	long double get_p_value(int v, hash_set<int>& sub_graph, hash_set<int>& around_sub_graph);
	long double get_C(int n, int k);
	void wirte2file();
	hash_set<int> get_around_node(hash_set<int>& sub_graph);
};