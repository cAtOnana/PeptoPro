#include"pFind_PairResearch.h"
#include<algorithm>
#include<unordered_map>
#include<fstream>
#include<iostream>
ostream & operator<<(ostream & os, const spectra & s)
{
	os << s.file_name << "	" << s.scan_no << "	" << s.exp_mh << "	" << s.charge << "	" << s.q_value << "	" << s.seq << "	"
		<< s.calc_mh << "	" << s.mass_shift << "	" << s.raw_score << "	" << s.final_score << "	" << s.modi << "	" << s.spec
		<< "	" << s.prot << "	" << s.posi << "	" << s.label << "	" << s.targe << "	" << s.mc_sites << "	" << s.afm_shift << "	"
		<< s.others;
	return os;
}

istream & operator>>(istream & in, vector<spectra> & list_result)
{
	char ch;
	spectra temp;//虽然是每次初始化，但在优化下会变成每次初始化都占用同一内存区域，从而使is_mut,modi等变量
	//直到赋值前都会保留上一次的值（因为内存恰好对上了）；因此对后续步骤有决定性作用的is_mut值，更应该用完整的if else语句重新赋值
	while (in >> temp.file_name) {
		in >> temp.scan_no;
		in >> temp.exp_mh;
		in >> temp.charge;
		in >> temp.q_value;
		in >> temp.seq;
		in >> temp.calc_mh;
		in >> temp.mass_shift;
		in >> temp.raw_score;
		in >> temp.final_score;
		in.get(ch).get(ch);
		if (ch == '\t')
			temp.modi = "";
		else
		{
			in >> temp.modi;
			temp.modi = ch + temp.modi;
			if (temp.modi.find(">") != string::npos) {
				temp.is_mut = true;
				//std::cout << temp.modi.find<< std::endl;
			}
			else
				temp.is_mut = false;
		}
		in >> temp.spec;
		in >> temp.prot;
		temp.prot=temp.prot.substr( temp.prot.find("ENSP"), 15);
		in >> temp.posi;
		in >> temp.label;
		in >> temp.targe;
		in >> temp.mc_sites;
		in >> temp.afm_shift;
		in >> temp.others;
		in >> temp.mut_count;
		list_result.push_back(temp);
	}
	return in;
}

int mark(vector<spectra>& list)//将数据分成几个大组，每组有相同的最小子序列，且不重复，但此函数未解决同时包含多条不同的最小子序列（从而可以归属于多个不同组）的问题。但愿数据集中没有这样的肽段序列
{
	sort(list.begin(), list.end(), sortbyleg);
	int mark = 1;//注意，mark从1开始。
	for (int i = 0; i < list.size(); i++) {
		if (list[i].marker != 0)
			continue;
		else
			list[i].marker = mark++;//赋值后更新
		for (int j = i + 1; j < list.size(); j++) {
			if (list[j].marker != 0)
				continue;
			if (list[j].seq.find(list[i].seq) != string::npos)//如果找到了的话
				list[j].marker = list[i].marker;
		}
	}
	sort(list.begin(), list.end(), sortbymarker);
	return mark;//返回mark，用于建立bool数组判断是否可输出。
}



bool sortbyleg(spectra & a, spectra & b)
{
	return a.seq.length()<b.seq.length();
}

bool sortbymarker(spectra & a, spectra & b)
{
	return a.marker<b.marker;
}

mut_pep_inform pepmutation(const spectra & p, ifstream& intri, unordered_map<string, char>& table)
{
	//文本处理
	string modi = p.modi;
	mut_pep_inform a;
		a.mutpep = p.seq;
		a.size = p.mut_count;
		a.pos_mut = new int[a.size]{ 0 }; 
	string section;
	while (modi.find(">") != string::npos && modi.find(";") != string::npos) {
		int end_sign = modi.find(";");
		section = modi.substr(0, end_sign + 1);
		modi.erase(0, end_sign + 1);
		for (int i = 0; section[i] != ','; i++)//得到突变坐标，存于pos_mut中
			*a.pos_mut = *a.pos_mut * 10 + section[i] - 49;
		section.erase(0, section.find(',') + 1);
		string ori_res, mut_res;
		ori_res=section.substr( 0, 3);
		int mut_sign = modi.find(">");
		mut_res=section.substr(mut_sign + 1, 3); 
		char ori, mut;
		if (table.find(ori_res) != table.end() && table.find(mut_res) != table.end()) {
			ori = table[ori_res];
			mut = table[mut_res];
		}
		else
			std::cout << "有找不到缩写的残基，其名为：" <<mut_res << std::endl;
		if (a.mutpep[*a.pos_mut] == ori)
			a.mutpep[*a.pos_mut] = mut;
		a.pos_mut++;
		ori_res.~basic_string();
		mut_res.~basic_string();
	}
	section.~basic_string();
	return a;
}

struct tri_hash{
	size_t operator()(const string & str)const {
		int _h = 0;
		for (int i = 0; i < str.length(); i++) {
			_h = _h * 2 + str[i] - 65;
		}
		return _h;
	}
};

struct tri_compare{
	bool operator()(const string& a1, const string& a2)const {
		return a1 == a2;
	}
};