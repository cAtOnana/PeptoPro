#include"reflect.h"
#include"pFind_PairResearch.h"
#include<vector>
#include<algorithm>
//result文件中肽段是以nm号标记的，非同义突变文件中的蛋白质序列是以ensp号标记的，现通过清理空白项后的对照表建立映射将两者对应起来后进行比对；无突变的肽段直接比对即可，有突变的肽段先在对应蛋白质序列中进行突变转换然后比对
using namespace std;
int argc = 5;
char* argv[5] = { "aaa","result.txt","非同义突变.txt","ENSP-NM.txt","Homo_sapiens.GRCh38.pep.all_cnew.fa" };
int main(){//int argc, char* argv[]) {//result mrna非同义突变 对照表 蛋白质序列fasta
	char* logname = "PeptoPro.log";
	ofstream log(logname);
	if (argc != 5) {
		cout << "参数数目异常，程序退出。\n";
		log << "参数数目异常，程序退出。\n";
		exit(EXIT_FAILURE);
	}
	ifstream inres(argv[1]);
	if (!inres.is_open()) {
		cout << argv[1] << "打开失败，程序退出\n";
		log << argv[1]<< "打开失败，程序退出\n";
		exit(EXIT_FAILURE);
	}
	ifstream inmrn(argv[2]);
	if (!inmrn.is_open()) {
		cout << argv[2] << "打开失败，程序退出\n";
		log << argv[2]<< "打开失败，程序退出\n";
		exit(EXIT_FAILURE);
	}
	ifstream incom(argv[3]);
	if (!incom.is_open()) {
		cout << argv[3] << "打开失败，程序退出\n";
		log << argv[3] << "打开失败，程序退出\n";
		exit(EXIT_FAILURE);
	}
	ifstream inref(argv[4]);
	if (!inref.is_open()) {
		cout << argv[4] << "打开失败，程序退出\n";
		log << argv[4]<< "打开失败，程序退出\n";
		exit(EXIT_FAILURE);
	}
	vector<spectra> list_result;
	inres >> list_result;
	int markcount =mark(list_result);
	for (int i = 0; i < list_result.size();i++)
		log << list_result[i];

	vector<pro> list_pro;
	inmrn >> list_pro;
	fillhseq(inref, list_pro);
	//输入对照表
	fillnm(incom, list_pro);
	//开始比对
	bool* cout_modi = new bool[markcount] {false};
	bool* cout_no = new bool[markcount] {false};
	///此两数组用于记录某一mark号对应的数据组能否输出，只有两者皆为true（即某一组中至少有一对经过验证的突变对）才能输出。故比对时也应按照此原则对应地更新数组。
	int saving = 0;
	///记录匹配到的list_pro的坐标，方便下一次匹配
	for (int i = 0; i < list_pro.size();i++) {
		if (list_result[0].is_modi == false) {
			if (list_pro[i].hseq.find(list_result[0].seq) != string::npos) {
				saving = i;
				cout_modi[list_result[0].marker] = true;
			}
		}
		else if (list_result[0].is_modi == true) {
			string temp = list_pro[i].hseq;
			temp.replace(list_pro[i].pos, 1,1 ,list_pro[i].mutataa);
			if (temp.find(list_result[0].seq ) != string::npos){
				saving = i;
				cout_no[list_result[0].marker] = true;
			}
		}
		
	}
	for (int i = 0; i < list_result.size();) {//更新放在分支后
		for (int j = 0; j < list_pro.size(); j++) {

		}
	}
		
}


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
		list_result.push_back(temp);
	}
	return in;
}

int mark(vector<spectra>& list)
{
	sort(list.begin(),list.end(), sortbyleg);
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
