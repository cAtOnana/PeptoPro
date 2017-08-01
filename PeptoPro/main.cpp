#include"reflect.h"
#include"pFind_PairResearch.h"
#include<vector>
#include<unordered_map>
//result文件中肽段是以nm号标记的，非同义突变文件中的蛋白质序列是以ensp号标记的，现通过清理空白项后的对照表建立映射将两者对应起来后进行比对；无突变的肽段直接比对即可，有突变的肽段先在对应蛋白质序列中进行突变转换然后比对
using namespace std;
int argc = 6;
char* argv[6] = { "aaa","result.txt","非同义突变.txt","ENSP-NM.txt","Homo_sapiens.GRCh38.pep.all_cnew.fa","三字表.txt" };
struct pro_hash {
	size_t operator()(const string& str) const {
		int _hash_value = 0;
		for (int i = 0; i < str.length(); i++)
			_hash_value = _hash_value * 2 + str[i] - 49;
		return _hash_value;
	}
};
struct pro_compare {
	bool operator()(const string& a1, const string& a2) const {
		return a1==a2;//试试此处不是相等判断而是包含判断时是否可行
	}
};
typedef unordered_map<string, int> ProIndexMap;
int main(){//int argc, char* argv[]) {//result mrna非同义突变 对照表 蛋白质序列fasta
	char* logname = "PeptoPro.log";
	ofstream log(logname);
	if (argc != 6) {
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
	//三字表 单子表文件名；用tab分隔
	ifstream intri(argv[5]);//输入三字表和单子表对应
	if (!intri.is_open()) {
		cout << argv[5] << "打开失败，程序退出\n";
		log << argv[5] << "打开失败，程序退出\n";
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
	ProIndexMap pimap;//NM号为key，index为value
	for (int i = 0; i < list_pro.size(); i++) {//填表
		//由于NM号不唯一且相互间以逗号间隔以双引号包围，故应先对其进行处理，剥出一个个单独的NM号，以形成一对一的 <NM号> = 序号 组合
		string temp = list_pro[i].nm;
		for (int k = 0; k < temp.length(); k++)
			if (temp[k] == '"')
				temp.erase(k, 1);
		while (temp.find(',')) {
			int g = temp.find(',');
			string str;
			str.assign(temp, 0, g);
			pimap[str] = i;
			temp.erase(0, g + 1);
		}
		pimap[temp] = i;//将最后一个NM号也录入表中
	}
	for (int i = 0; i < list_result.size();) {//更新放在分支末
		spectra& pep = list_result[i];//起个简单的别名方便编程，提高可读性
		if (pep.marker != mark_saving) {
			if (pep.is_mut == false) {
				if (pimap.find(pep.prot) != pimap.end()) {//找到了的话
					int j = pimap[pep.prot];
					if (list_pro[j].nm.find(pep.prot) != string::npos) {
						pos_saving = j;//传递给pos_saving
						mark_saving =pep.marker;
						cout_no[mark_saving] = true;
						i++;//update
						break;
					}
				}
					else {//找不到的话，擦除，然后继续循环
						vector<spectra>::iterator itor = list_result.begin() + i;
						list_result.erase(itor);//erase & update
					}
			}
			else {//类似上面
				if (pimap.find(pep.prot)!=pimap.end()) {//先找后改，高效
					int j = pimap[pep.prot];
					string temp = list_pro[j].hseq;
					temp = temp.replace(list_pro[j].pos, 1, 1, list_pro[j].mutataa);
					mut_pep_inform temp2 = pepmutation(pep,intri);//突变肽链
					auto pos_find = temp.find(temp2.mutpep);
					bool access = false;//此变量用于判定是否有至少一个肽段突变位点与蛋白质突变位点重合
					for (int i = 0; i < temp2.size; i++) {
						if (pos_find + temp2.pos_mut[i] == list_pro[j].pos)
							access = true;
					}
					if (pos_find != string::npos && access) {//此处pos_mut是result中肽段突变位点的坐标，与匹配起始位点坐标（pos_find）相加应等于蛋白质序列上突变位点坐标
						pos_saving = j;
						mark_saving = pep.marker;
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
		else {//mark_saving与pep.marker相同时
			if (pep.is_mut == false) {
				if (list_pro[pos_saving].hseq.find(pep.seq) != string::npos) {
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
				if (list_pro[pos_saving].nm.find(pep.prot) != string::npos) {
					string temp = list_pro[pos_saving].hseq;
					temp = temp.replace(list_pro[pos_saving].pos, 1, 1, list_pro[pos_saving].mutataa);//突变蛋白链
					mut_pep_inform temp2 = pepmutation(pep,intri);//突变肽链
					auto pos_find = temp.find(temp2.mutpep);
					bool access = false;//此变量用于判定是否有至少一个肽段突变位点与蛋白质突变位点重合
					for (int i = 0; i < temp2.size; i++) {
						if (pos_find + temp2.pos_mut[i] == list_pro[pos_saving].pos)
							access = true;
					}
					if (pos_find != string::npos && access) {//此处pos_mut是result中肽段突变位点的坐标，与匹配起始位点坐标（pos_find）相加应等于蛋白质序列上突变位点坐标
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
	/////输出
	string outname = "OpenResearch筛选后.txt";
	ofstream out(outname);
	for (int i = 0; i < list_result.size(); i++) {
		if (cout_modi[list_result[i].marker] && cout_no[list_result[i].marker])
			out << list_result[i];
	}
	cout << "Done!" << endl;
	char a;
	cin >> a;
	return 0;
}

