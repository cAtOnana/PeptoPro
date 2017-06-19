#include "reflect.h"

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
