#include "graph_cluster.h"

int GraphCluster::max_degree_unvisited()
{
	int value = -1;
	int index = -1;
	for (int i = 0; i < n; i++)
	{
		if (!visited[i] && ki[i]>value)
		{
			value = ki[i];
			index = i;
		}
	}
	return index;
}

hash_set<int> GraphCluster::get_around_node(hash_set<int>& sub_graph)
{
	DS = 0;
	hash_set<int> around_sub_graph;
	for (hash_set<int>::iterator it = sub_graph.begin(); it != sub_graph.end(); it++)
	{
		int u = *it;
		DS += ki[u];
		for (int j = 0; j < graph[u].size(); j++)
		{
			int v = graph[u][j];
			if (sub_graph.count(v) == 0 && around_sub_graph.count(v) == 0)
			{
				around_sub_graph.insert(v);
			}
		}
	}
	return around_sub_graph;
}

void GraphCluster::SSF_old()
{
	int v = max_degree_unvisited();
	int cnt = 0;
	while (v != -1)
	{
		visited[v] = true;
		sub_graph_node = graph[v];
		sub_graph_node.push_back(v);
		vector<int> temp_S = search_one_complex();

		if (temp_S.size() > 2)
		{
			S.push_back(temp_S);
			
			for (int i = 0; i < temp_S.size(); i++)
			{
				if (!visited[temp_S[i]])
					cnt++;
				visited[ temp_S[i] ] = true;
			}
			
		}
		v = max_degree_unvisited();

		cnt += 1;
		cout << cnt << "  " << temp_S.size()<< endl;

	}
	
	Setcover();
	wirte2file();
}

void GraphCluster::SSF()
{
	
	
	//Setcover();
	wirte2file();
}

vector<int> GraphCluster::search_one_complex()
{
	hash_set<int>sub_graph;
	hash_set<int>around_sub_graph;
	for (int i = 0; i < sub_graph_node.size(); i++)
	{
		int u = sub_graph_node[i];
		sub_graph.insert(u);
	}

	around_sub_graph = get_around_node(sub_graph);

	hash_set<int> last_sub_graph = sub_graph;
	int cnt = 0;
	while (true)
	{
		cnt++;
		hash_map<int, long double> p_value;
		for (hash_set<int>::iterator it = sub_graph.begin(); it != sub_graph.end(); it++)
		{
			int v = *it;
			p_value[v] = get_p_value(v, sub_graph, around_sub_graph);
		}
		for (hash_set<int>::iterator it = around_sub_graph.begin(); it != around_sub_graph.end(); it++)
		{
			int v = *it;
			p_value[v] = get_p_value(v, sub_graph, around_sub_graph);
		}

		hash_set<int> new_sub_graph;
		for (hash_map<int, long double>::iterator it = p_value.begin(); it != p_value.end(); it++)
		{
			if (it->second <= alpha / n)
			{
				new_sub_graph.insert(it->first);
			}
		}
		
		if (new_sub_graph.size() == 0 || cnt >= 30 || is_same_graph(new_sub_graph, sub_graph))
		{
			sub_graph = new_sub_graph;
			break;
		}
		sub_graph = new_sub_graph;
		around_sub_graph = get_around_node(sub_graph);	
	}

	vector<int> ret;
	for (hash_set<int>::iterator it = sub_graph.begin(); it != sub_graph.end(); it++)
	{
		ret.push_back(*it);
	}
	return ret;
}

bool GraphCluster::is_same_graph(hash_set<int> &graph1, hash_set<int>&graph2)
{
	if (graph1.size() != graph2.size())
	{
		return false;
	}
	for (hash_set<int>::iterator it = graph1.begin(); it != graph1.end(); it++)
	{
		if (graph2.count(*it) == 0)
		{
			return false;
		}
	}
	return true;
}

void GraphCluster::Setcover()
{
	F.clear();
	hash_set<int> W;
	vector<bool> used(S.size(), false);
	for (int i = 0; i < S.size(); i++)
	{
		for (int j = 0; j < S[i].size(); j++)
		{
			W.insert(S[i][j]);
		}
	}
	while (W.size() > 0)
	{
		int value = -1;
		int index = -1;
		for (int i = 0; i < S.size(); i++)
		{
			if (used[i]) continue;
			int temp = 0;
			for (int j = 0; j < S[i].size(); j++)
			{
				if (W.count(S[i][j]) >0)
				{
					temp++;
				}
			}
			if (temp > value)
			{
				value = temp;
				index = i;
			}
		}
		if (index != -1)
		{
			used[index] = true;
		}
		else
		{
			break;
		}
		F.push_back(S[index]);
		for (int i = 0; i < S[index].size(); i++)
		{
			W.erase(S[index][i]);
		}
	}
}