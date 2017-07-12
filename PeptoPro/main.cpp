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
	int pos_saving = 0, mark_saving = 0;
	///pos_saving记录匹配到的list_pro的坐标，方便下一次匹配；mark_saving用于记录当前组号
	///mark_saving记录当前正在执行匹配的肽序列所属的组；同时用于判定pos_saving的使用与否，解决边界问题（因为pos_saving是进入直接比对流程的）
	for (int i = 0; i < list_result.size();) {//更新放在分支末
		if (list_result[i].marker != mark_saving) {
			if (list_result[i].is_modi == false) {
				int j = 0;
				for (; j < list_pro.size(); j++) {
					if (list_pro[j].nm.find(list_result[i].prot)!=string::npos && list_pro[j].hseq.find(list_result[i].seq) != string::npos) {//找到了的话
						pos_saving = j;
						mark_saving = list_result[i].marker;
						cout_no[mark_saving] = true;
						i++;//update
						break;
					}
					else {//找不到的话，擦除，然后继续循环
						vector<spectra>::iterator itor = list_result.begin() + i;
						list_result.erase(itor);//erase & update
					}
				}
			}
			else {//类似上面
				int j = 0;
				for (; j < list_pro.size(); j++) {
					string temp = list_pro[j].hseq;
					temp = temp.replace(list_pro[j].pos, 1, 1, list_pro[j].mutataa);//A ?在目标肽链上有多个突变，怎么办？
					string temp2 = list_result[i].prot;
					//此处将temp2改为突变序列
					if (list_pro[j].nm.find(list_result[i].prot) != string::npos) {
						auto pos_find = temp.find(temp2);
						if (pos_find != string::npos && pos_find + pos_mut == list_pro[j].pos) {//此处pos_mut是result中肽段突变位点的坐标，与匹配起始位点坐标（pos_find）相加应等于蛋白质序列上突变位点坐标
							pos_saving = j;
							mark_saving = list_result[i].marker;
							cout_modi[mark_saving] = true;
							i++;//update
							break;
						}
					}
					else {
						vector<spectra>::iterator itor = list_result.begin() + i;
						list_result.erase(itor);//erase & update
					}
				}
			}
		}
		else {
			if (list_result[i].is_modi == false) {
				if (list_pro[pos_saving].hseq.find(list_result[i].seq) != string::npos) {
					cout_no[mark_saving] = true;
					i++;//update
					break;
				}
				else {
					vector<spectra>::iterator itor = list_result.begin()+i;
					list_result.erase(itor);//erase & update
				}
			}
			else {
				string temp = list_pro[pos_saving].hseq;
				temp = temp.replace(list_pro[pos_saving].pos, 1, 1, list_pro[pos_saving].mutataa);//问题A
				string temp2 = list_result[i].prot;
				//此处将temp2改为突变序列
				if (list_pro[pos_saving].nm.find(list_result[i].prot) != string::npos) {
					auto pos_find = temp.find(temp2);
					if (pos_find != string::npos && pos_find + pos_mut == list_pro[j].pos) {//此处pos_mut是result中肽段突变位点的坐标，与匹配起始位点坐标（pos_find）相加应等于蛋白质序列上突变位点坐标
						cout_modi[mark_saving] = true;
						i++;//update
						break;
					}
				}
				else {
					vector<spectra>::iterator itor = list_result.begin() + i;
					list_result.erase(itor);// erase & update
				}
			}
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

int mark(vector<spectra>& list)//将数据分成几个大组，每组有相同的最小子序列，且不重复，但此函数未解决同时包含多条不同的最小子序列（从而可以归属于多个不同组）的问题。但愿数据集中没有这样的肽段序列
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
