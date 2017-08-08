#include"reflect.h"
#include"pFind_PairResearch.h"
#include<vector>
#include<unordered_map>
#include<algorithm>
//result文件中肽段是以nm号标记的，非同义突变文件中的蛋白质序列是以ensp号标记的，现通过清理空白项后的对照表建立映射将两者对应起来后进行比对；无突变的肽段直接比对即可，有突变的肽段先在对应蛋白质序列中进行突变转换然后比对
using namespace std;
int argc = 6;
char* argv[6] = { "aaa","all_openresearch_data_esemble_pairfinding_result.txt","非同义突变.txt","ENSP-NM.txt","Homo_sapiens.GRCh38.pep.all_cnew.fa","三字表.txt" };
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
typedef unordered_map<string, vector<int>> ProIndexMap;
bool comparestr(pro& a, pro& b) {
	return a.ensp < b.ensp;
}
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

	vector<pro> list_pro;
	inmrn >> list_pro;
	fillhseq(inref, list_pro);
	//开始比对
	///pos_saving记录匹配到的list_pro的坐标，方便下一次匹配；mark_saving用于记录当前组号
	///mark_saving记录当前正在执行匹配的肽序列所属的组；同时用于判定pos_saving的使用与否，解决边界问题（因为pos_saving是进入直接比对流程的）
	ProIndexMap pimap;//ENSP号为key，index为value
	sort(list_pro.begin(), list_pro.end(), comparestr);
	string ensp = list_pro[0].ensp;
	vector<int> index;
	for (int i = 0; i < list_pro.size(); i++) {//填表
		index.push_back(i);
		if (i == list_pro.size() - 1) {
			pimap[ensp] = index;
		}
		if (list_pro[i].ensp != list_pro[i + 1].ensp) {
			pimap[ensp] = index;
			ensp = list_pro[i + 1].ensp;
			index.clear();
		}
	}
	//构建三字表-单字表映射
	unordered_map<string, char> table;
	string tri;
	char mono;
	while (intri >> tri) {
		intri >> mono;
		table[tri] = mono;
	}
	bool* cout_modi = new bool[markcount] {false};
	bool* cout_no = new bool[markcount] {false};
	///此两数组用于记录某一mark号对应的数据组能否输出，只有两者皆为true（即某一组中至少有一对经过验证的突变对）才能输出。故比对时也应按照此原则对应地更新数组。
	int mark_saving = 0;
	for (int i = 0; i < list_result.size();i++) {
		spectra& pep = list_result[i];//起个简单的别名方便编程，提高可读性
		//if (pep.marker != mark_saving) {
			if (pep.is_mut == false) {
				if (pimap.find(pep.prot) != pimap.end()) {//找到了的话
					vector<int> j = pimap[pep.prot];
					//pos_saving = j;//传递给pos_saving
					mark_saving =pep.marker;
					for (int i = 0; i < j.size(); i++) {
						if (list_pro[j[i]].hseq.find(pep.seq) != string::npos) {
							cout_no[mark_saving] = true;
							break;
						}
						else
							list_result[i].outputable = false;
					}					
				}
					else {
						list_result[i].outputable = false;
						//找不到的话，擦除，然后继续循环
					}
			}
			else {//类似上面
				if (pimap.find(pep.prot)!=pimap.end()) {//先找后改，高效
					vector<int> j = pimap[pep.prot];
					//当处理openresearch数据时，启用暗绿语句
					///mut_pep_inform temp2 = pepmutation(pep,intri,table);//突变肽链
					for (int i = 0; i < j.size(); i++) {//循环突变蛋白质序列并比对
						string temp = list_pro[j[i]].hseq;
						temp = temp.replace(list_pro[j[i]].pos, 1, 1, list_pro[j[i]].mutataa);
						///auto pos_find = temp.find(temp2.mutpep);
					//当处理理论突变数据时，启用暗绿语句
						auto pos_find = temp.find(pep.seq);
						///bool access = false;//此变量用于判定是否有至少一个肽段突变位点与蛋白质突变位点重合
						///for (int i = 0; i < temp2.size; i++) {//此处pos_mut是result中肽段突变位点的坐标(0起始)，与匹配起始位点坐标（pos_find）（同样是0起始）相加应等于蛋白质序列上突变位点坐标（0起始）
						///	if (pos_find + temp2.pos_mut[i] == list_pro[j[i]].pos)
						///		access = true;
						///}
						if (pos_find != string::npos){ ///&& access) {
							mark_saving = pep.marker;
							cout_modi[mark_saving] = true;
							break;
						}
						else 
							list_result[i].outputable = false;
					}
				}
				else {
					list_result[i].outputable = false;
				}
			}
		}
	/////输出
	string outname = "openresearch筛选后.txt";
	ofstream out(outname);
	for (int i = 0; i < list_result.size(); i++) {
		if (cout_modi[list_result[i].marker] && cout_no[list_result[i].marker]&& list_result[i].outputable)
			out << list_result[i]<<endl;
	}
	cout << "Done!" << endl;
	char a;
	cin >> a;
	return 0;
}

