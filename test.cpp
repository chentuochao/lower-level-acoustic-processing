#include <vector>
#include <iostream>
#include<string.h>
using namespace std;

int main( )
{

vector<int> v1, v2, v3;
vector<int>::iterator iter;
v1.reserve(100);
v2.reserve(100);

v1.push_back(10);
v1.push_back(20);
v1.push_back(30);
v1.push_back(40);
v1.push_back(50);
v2.push_back(1);
v2.push_back(2);

cout << "v1 = " ;
for (iter = v1.begin(); iter != v1.end(); iter++)
cout << *iter << " ";
cout << endl;

cout << "v2 = "<< v2.capacity() << ' ';
for (iter = v2.begin(); iter != v2.end(); iter++)
cout << *iter << " ";
cout << endl;
cout << "ptr address: " << (void*)v2.data() << endl;

/*
v2 = v1;
cout << "v2 = ";
for (iter = v2.begin(); iter != v2.end(); iter++)
cout << *iter << " ";
cout << endl;
cout << "ptr address: " << (void*)v2.data() << endl;
*/

v2.assign(v1.begin(), v1.end());
cout << "v2 = " << v2.capacity() << ' ';
for (iter = v2.begin(); iter != v2.end(); iter++)
cout << *iter << " ";
cout << endl;
cout << "ptr address: " << (void*)v2.data() << endl;

int *p2 = v2.data();


int * p = new int[v2.size()];
//copy(v2.begin(), v2.end(), p);
memcpy(p, p2, sizeof(int)*v2.size());

vector<int>::iterator it1, it2;
it1 = v2.begin();
it2 = v2.end() - 2;
v2.erase(it1, it2);
cout << "v2 = " << v2.capacity() << ' ';
for (iter = v2.begin(); iter != v2.end(); iter++)
cout << *iter << " ";
cout << endl;
cout << "ptr address: " << (void*)v2.data() << endl;

cout << "p = ";
for (int i = 0; i < v2.size(); i++)
cout << p[i] << " ";
cout << endl;

cout << "v1 address: " << (void*)v1.data() << endl;
v1.clear() ;
v1.push_back(0);
cout << "v1 = ";
for (iter = v1.begin(); iter != v1.end(); iter++)
cout << *iter << " ";
cout << endl;
cout << "v1 address: " << (void*)v1.data() << endl;

return 0;
}