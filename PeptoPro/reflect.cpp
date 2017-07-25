#include "reflect.h"
#include<algorithm>
#include<unordered_map>
struct str_hash {//自定义hash函数
	size_t operator() (const string& str) const {
		unsigned long _h = 0;
		for (int i = 0; i < str.size(); i++)
			_h = _h * 2 + str[i] - 49;
		return size_t(_h);
	}
};
struct str_compare {
	bool operator()(const string& a1, const string& a2) const {
		return a1 == a2;
	}
};
typedef unordered_map<string, string, str_hash, str_compare> NmEnspMap,EnspHseqMap;
//////////////////////////////////////////////////////////////
istream & operator>>(istream & is, vector<pro> list)
{
	pro temp;
	string waste;
	getline(is, waste);
	getline(is, waste);
	while (is >> waste) {
		for (int i = 0; i < 12; i++)
			is >> waste;
		is >> temp.ensp;
		is >> temp.origaa;
		is >> temp.mutataa;
		is >> temp.pos;
	}
	return is;
}
bool sortbyensp(pro & a, pro & b)
{
	return a.ensp < b.ensp;
}

void fillhseq(istream & is, vector<pro> list)
{
	//sort(list.begin(), list.end(), sortbyensp);
	//string temp;
	//int a = 0;
	//int b = list.size();
	//int mid = (a + b) / 2;
	//while (is.get()=='>') {
	//	is >> temp;
	//	while (a < b) {//二分法迅速确认位置
	//		if (list[mid].ensp == temp) {
	//			is >> list[mid].hseq;
	//			break;
	//		}
	//		else if (temp > list[mid].ensp) {
	//			a = mid;
	//			mid = (a + b) / 2;
	//		}
	//		else if (temp < list[mid].ensp) {
	//			b = mid;
	//			mid = (a + b) / 2;
	//		}
	//	}
	//}
	////以上为老方法，下面尝试hash_map法
	EnspHseqMap ehmap;
	string ensp, hseq;
	while (is) {//填表
		while (is.get() != '>')//找'>'
			continue;
		is >> ensp >> hseq;
		ehmap[ensp] = hseq;
	}
	for (int i = 0; i < list.size(); i++) {
		if (ehmap.find(list[i].ensp) != ehmap.end())
			list[i].hseq = ehmap[list[i].ensp];
	}
}

void fillnm(istream & is, vector<pro> list)
{
	//sort(list.begin(), list.end(), sortbyensp);
	//string temp;
	//int a = 0;
	//int b = list.size();
	//int mid = (a + b) / 2;
	//while (is >> temp) {
	//	while (a<b) {//二分法迅速确认位置
	//		if (list[mid].ensp == temp) {
	//			is >> list[mid].nm;
	//			break;
	//		}
	//		else if (temp>list[mid].ensp) {
	//			a = mid;
	//			mid = (a + b) / 2;
	//		}
	//		else if (temp < list[mid].ensp) {
	//			b = mid;
	//			mid = (a + b) / 2;
	//		}
	//	}
	/////以上为老方法，下面尝试hash_map法

	string nmname, enspname;
	NmEnspMap mymap;
	string nm_name, ensp_name;
	while (is >> ensp_name) {//填表
		is >> nm_name;
		mymap[ensp_name] = nm_name;
	}
	for (int i = 0; i < list.size(); i++){
		if (mymap.find(list[i].ensp) != mymap.end())
			list[i].nm = mymap[list[i].ensp];
		else
			cout << "没有找到" << list[i].ensp << "相对应的NM号。";
	}
}
