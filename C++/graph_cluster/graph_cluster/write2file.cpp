#include "graph_cluster.h"

void GraphCluster::wirte2file()
{
	ofstream fout(out_file);
	cout << out_file << endl;
	for (int i = 0; i < F.size(); i++)
	{
		fout << int2string[F[i][0]];
		for (int j = 1; j < F[i].size(); j++)
		{
			fout << " " << int2string[F[i][j]];
		}
		fout << endl;
	}
	fout.close();
}