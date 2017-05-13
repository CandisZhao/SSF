#include "graph_cluster.h"

long double GraphCluster::get_C(int n, int k)
{
	return f[n] - f[k] - f[n - k];
}

long double GraphCluster::get_p_value(int u, hash_set<int>& sub_graph, hash_set<int>& around_sub_graph)
{
	if (sub_graph.size() == 1)
	{
		return 1;
	}

	int C[2][2] = {0};
	for (int i = 0; i < graph[u].size(); i++)
	{
		int v = graph[u][i];
		if (sub_graph.count(v) != 0)
		{
			C[0][0] += 1;
		}
	}
	C[1][0] = ki[u] - C[0][0];
	if (sub_graph.count(u) != 0)
	{
		C[0][1] = DS - ki[u] - C[0][0];
		C[1][1] = D - DS - C[1][0];
	}
	else
	{
		C[0][1] = DS - C[0][0];
		C[1][1] = D - DS - ki[u] - C[1][0];
	}
	long double P = get_C(C[0][0] + C[0][1], C[0][0]) + get_C(C[1][0] + C[1][1], C[1][0]) - get_C(D - ki[u], ki[u]);
	P = exp(P);
	P = P > 1.0? 1.0: P;

	long double q1 = 1.0*C[1][0] * C[0][1] / (C[0][0]*C[1][1]);

	if (q1 >= 1- 1e-10)
	{
		return 1.0;
	}
	else
	{
		long double rou = C[0][0] * (D - ki[u]) / (ki[u] * (C[0][0] + C[0][1]));
		if (rou <= 2)
		{
			return P / (1 - q1);
		}
		else
		{
			long double value = P / (1.0 - q1) * (1-pow(q1, C[1][0]+1) );
			return value;
		}
	}
}