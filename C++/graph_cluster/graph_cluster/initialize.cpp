#include "graph_cluster.h"

GraphCluster::GraphCluster(long double alp, string input, string output)
{
	n = 0;
	D = 0;
	alpha = alp;
	in_file = input;
	out_file = output;
}

void GraphCluster::initialize()
{
	ifstream fin(in_file);
	vector<string> lines;
	vector<string> pools;
	string str;
	vector<int> v_t;
	while (getline(fin, str))
	{
		lines.push_back(str);
		stringstream ss(str);
		string  temp[2];
		ss >> temp[0] >> temp[1];
		for (int i = 0; i < 2; i++)
		{
			if (string2int.count(temp[i]) == 0)
			{
				string2int[temp[i]] = -1;
				pools.push_back(temp[i]);
				graph.push_back(v_t);
			}
		}
	}
	fin.close();
	
	sort(pools.begin(), pools.end());
	for (int i = 0; i < pools.size(); i++)
	{
		string2int[pools[i]] = n;
		int2string[n++] = pools[i];
	}

	for (int i = 0; i < lines.size(); i++)
	{
		string str = lines[i];
		stringstream ss(str);
		string  temp[2];
		ss >> temp[0] >> temp[1];

		int u = string2int[temp[0]];
		int v = string2int[temp[1]];
		graph[u].push_back(v);
		graph[v].push_back(u);
	}

	vector<bool> temp(n, false);
	visited = temp;
	for (int i = 0; i < n; i++)
	{
		sort(graph[i].begin(), graph[i].end());
		graph[i].erase(unique(graph[i].begin(), graph[i].end()), graph[i].end());

		ki.push_back(graph[i].size());
		D += graph[i].size();
	}

	long double value = 0;
	f.push_back(value);
	for (int i = 1; i <= D; i++)
	{
		value = value + log(1.0*i);
		f.push_back(value);
	}
}

void GraphCluster::initialize2()
{
	ifstream fin(in_file);
	string str;
	vector<int> v_t;
	while (getline(fin, str))
	{
		stringstream ss(str);
		string  temp[2];
		ss >> temp[0] >> temp[1];
		for (int i = 0; i < 2; i++)
		{
			if (string2int.count(temp[i]) == 0)
			{
				string2int[temp[i]] = n;
				int2string[n++] = temp[i];
				graph.push_back(v_t);
			}
		}
		int u = string2int[temp[0]];
		int v = string2int[temp[1]];
		graph[u].push_back(v);
		graph[v].push_back(u);
	}
	fin.close();
	
	vector<bool> temp(n, false);
	visited = temp;
	for (int i = 0; i < n; i++)
	{
		sort(graph[i].begin(), graph[i].end());
		graph[i].erase(unique(graph[i].begin(), graph[i].end()), graph[i].end());

		ki.push_back(graph[i].size());
		D += graph[i].size();
	}

	long double value = 0;
	f.push_back(value);
	for (int i = 1; i <= D; i++)
	{
		value = value + log(1.0*i);
		f.push_back(value);
	}
}