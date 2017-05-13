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
	
	long double alpha = 0.01;    // ����
	string in_file;
	string out_file;

	vector<long double> f;


	int n;                        // ͼ�ڵ���Ŀ
	int D;
	int DS;
	hash_map<string, int> string2int;
	hash_map<int, string> int2string;

	vector< vector <int> > graph;   // �ڽӱ�
	vector<int> sub_graph_node;
	vector<hash_set<int>> clusters;     // ��ͼ�ڵ�

	vector<int> ki;                 // �ڵ��
	vector<int> ki_in;              // �ڵ�����ͼ���ڶ�
	vector<int> ki_out;             // �ڵ��-�ڵ�����ͼ���ڶ�

	vector< vector <int> > S;
	vector< vector <int> > F;
	vector<bool> visited;           // �ѱ����Ľڵ�


	GraphCluster(long double alp, string input, string output);
	void initialize();              // ��ʼ��
	void initialize2();

	void SSF();   // ������
	void SSF_old();   // ������
	vector<int> search_one_complex();  
	int max_degree_unvisited();    
	void Setcover();                // ����

	bool is_same_graph(hash_set<int> &graph1, hash_set<int>&graph2);
	long double get_p_value(int v, hash_set<int>& sub_graph, hash_set<int>& around_sub_graph);
	long double get_C(int n, int k);
	void wirte2file();
	hash_set<int> get_around_node(hash_set<int>& sub_graph);
};