#include "conf.h"
#include "graph_cluster.h"
#include <cstdlib>
#include <ctime>
#include <set>
using namespace std;

int main()
{
	GraphCluster gc(alpha, input_file, output_file);
	gc.initialize();
	cout << gc.n << "  " << gc.D << endl;
	gc.SSF_old();
	return 0;
}