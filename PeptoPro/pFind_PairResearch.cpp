#include"pFind_PairResearch.h"
#include<algorithm>
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
	spectra temp;
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
			temp.is_modi = true;
		}
		in >> temp.spec;
		in >> temp.prot;
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

int mark(vector<spectra>& list)//�����ݷֳɼ������飬ÿ������ͬ����С�����У��Ҳ��ظ������˺���δ���ͬʱ����������ͬ����С�����У��Ӷ����Թ����ڶ����ͬ�飩�����⡣��Ը���ݼ���û���������Ķ�����
{
	sort(list.begin(), list.end(), sortbyleg);
	int mark = 1;//ע�⣬mark��1��ʼ��
	for (int i = 0; i < list.size(); i++) {
		if (list[i].marker != 0)
			continue;
		else
			list[i].marker = mark++;//��ֵ�����
		for (int j = i + 1; j < list.size(); j++) {
			if (list[j].marker != 0)
				continue;
			if (list[j].seq.find(list[i].seq) != string::npos)//����ҵ��˵Ļ�
				list[j].marker = list[i].marker;
		}
	}
	sort(list.begin(), list.end(), sortbymarker);
	return mark;//����mark�����ڽ���bool�����ж��Ƿ�������
}



bool sortbyleg(spectra & a, spectra & b)
{
	return a.seq.length()<b.seq.length();
}

bool sortbymarker(spectra & a, spectra & b)
{
	return a.marker<b.marker;
}

mut_pep_inform pepmutation(const spectra & p)
{
	//�ı�����
	string modi = p.modi;
	mut_pep_inform a;
		a.mutpep = p.seq;
		a.size = p.mut_count;
		a.pos_mut = new int[a.size]{ 0 }; 
	string section;
	while (modi.find("->") != string::npos && modi.find(';') != string::npos) {
		int end_sign = modi.find(';'), mut_sign = modi.find('>');
		section.assign(modi, 0, end_sign + 1);
		modi.erase(0, end_sign + 1);
		for (int i = 0; section[i] != ','; i++)//�õ�ͻ�����꣬����pos_mut��
			*a.pos_mut = *a.pos_mut * 10 + section[i] - 49;
		section.erase(0, section.find(',') + 1);
		string ori_res, mut_res;
		ori_res.assign(section, 0, 3);
		mut_res.assign(section, mut_sign + 1, 3);
		char ori, mut;
		////�ӹ�ϣ�����ҵ����ַ��л���Ӧ
		if (a.mutpep[*a.pos_mut] == ori)
			a.mutpep[*a.pos_mut] = mut;
		a.pos_mut++;
	}
	return a;
}
