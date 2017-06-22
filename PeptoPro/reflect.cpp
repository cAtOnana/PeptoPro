#include "reflect.h"
#include<algorithm>

istream & operator>>(istream & is, vector<pro> list)
{
	pro temp;
	string waste;
	getline(is, waste);
	getline(is, waste);
	while (is >> waste) {
		for(int i=0;i<12;i++)
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
	return a.ensp<b.ensp;
}

void fillhseq(istream & is, vector<pro> list)
{
	sort(list.begin(), list.end(), sortbyensp);
	string temp;
	int a = 0;
	int b = list.size();
	int mid = (a + b) / 2;
	while(is >> temp) {
		while (a<b) {//二分法迅速确认位置
			if (list[mid].ensp == temp) {
				is >> list[mid].hseq;
				break;
			}
			else if(temp>list[mid].ensp){
				a = mid;
				mid = (a + b) / 2;
			}
			else if (temp < list[mid].ensp) {
				b = mid;
				mid = (a + b) / 2;
			}
		}
	}
}

void fillnm(istream & is, vector<pro> list)
{
	sort(list.begin(), list.end(), sortbyensp);
	string temp;
	int a = 0;
	int b = list.size();
	int mid = (a + b) / 2;
	while (is >> temp) {
		while (a<b) {//二分法迅速确认位置
			if (list[mid].ensp == temp) {
				is >> list[mid].nm;
				break;
			}
			else if (temp>list[mid].ensp) {
				a = mid;
				mid = (a + b) / 2;
			}
			else if (temp < list[mid].ensp) {
				b = mid;
				mid = (a + b) / 2;
			}
		}
	}
}
