#include"reflect.h"
#include"pFind_PairResearch.h"
#include<vector>
#include<unordered_map>
#include<algorithm>
//result�ļ����Ķ�����nm�ű�ǵģ���ͬ��ͻ���ļ��еĵ�������������ensp�ű�ǵģ���ͨ������հ����Ķ��ձ���ӳ�佫���߶�Ӧ��������бȶԣ���ͻ����Ķ�ֱ�ӱȶԼ��ɣ���ͻ����Ķ����ڶ�Ӧ�����������н���ͻ��ת��Ȼ��ȶ�
using namespace std;
int argc = 6;
char* argv[6] = { "aaa","all_openresearch_data_esemble_pairfinding_result.txt","��ͬ��ͻ��.txt","ENSP-NM.txt","Homo_sapiens.GRCh38.pep.all_cnew.fa","���ֱ�.txt" };
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
		return a1==a2;//���Դ˴���������ж϶��ǰ����ж�ʱ�Ƿ����
	}
};
typedef unordered_map<string, vector<int>> ProIndexMap;
bool comparestr(pro& a, pro& b) {
	return a.ensp < b.ensp;
}
int main(){//int argc, char* argv[]) {//result mrna��ͬ��ͻ�� ���ձ� ����������fasta
	char* logname = "PeptoPro.log";
	ofstream log(logname);
	if (argc != 6) {
		cout << "������Ŀ�쳣�������˳���\n";
		log << "������Ŀ�쳣�������˳���\n";
		exit(EXIT_FAILURE);
	}
	ifstream inres(argv[1]);
	if (!inres.is_open()) {
		cout << argv[1] << "��ʧ�ܣ������˳�\n";
		log << argv[1]<< "��ʧ�ܣ������˳�\n";
		exit(EXIT_FAILURE);
	}
	ifstream inmrn(argv[2]);
	if (!inmrn.is_open()) {
		cout << argv[2] << "��ʧ�ܣ������˳�\n";
		log << argv[2]<< "��ʧ�ܣ������˳�\n";
		exit(EXIT_FAILURE);
	}
	ifstream incom(argv[3]);
	if (!incom.is_open()) {
		cout << argv[3] << "��ʧ�ܣ������˳�\n";
		log << argv[3] << "��ʧ�ܣ������˳�\n";
		exit(EXIT_FAILURE);
	}
	ifstream inref(argv[4]);
	if (!inref.is_open()) {
		cout << argv[4] << "��ʧ�ܣ������˳�\n";
		log << argv[4]<< "��ʧ�ܣ������˳�\n";
		exit(EXIT_FAILURE);
	}
	//���ֱ� ���ӱ��ļ�������tab�ָ�
	ifstream intri(argv[5]);//�������ֱ�͵��ӱ��Ӧ
	if (!intri.is_open()) {
		cout << argv[5] << "��ʧ�ܣ������˳�\n";
		log << argv[5] << "��ʧ�ܣ������˳�\n";
		exit(EXIT_FAILURE);
	}
	vector<spectra> list_result;
	inres >> list_result;
	int markcount =mark(list_result);

	vector<pro> list_pro;
	inmrn >> list_pro;
	fillhseq(inref, list_pro);
	//��ʼ�ȶ�
	///pos_saving��¼ƥ�䵽��list_pro�����꣬������һ��ƥ�䣻mark_saving���ڼ�¼��ǰ���
	///mark_saving��¼��ǰ����ִ��ƥ����������������飻ͬʱ�����ж�pos_saving��ʹ����񣬽���߽����⣨��Ϊpos_saving�ǽ���ֱ�ӱȶ����̵ģ�
	ProIndexMap pimap;//ENSP��Ϊkey��indexΪvalue
	sort(list_pro.begin(), list_pro.end(), comparestr);
	string ensp = list_pro[0].ensp;
	vector<int> index;
	for (int i = 0; i < list_pro.size(); i++) {//���
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
	//�������ֱ�-���ֱ�ӳ��
	unordered_map<string, char> table;
	string tri;
	char mono;
	while (intri >> tri) {
		intri >> mono;
		table[tri] = mono;
	}
	bool* cout_modi = new bool[markcount] {false};
	bool* cout_no = new bool[markcount] {false};
	///�����������ڼ�¼ĳһmark�Ŷ�Ӧ���������ܷ������ֻ�����߽�Ϊtrue����ĳһ����������һ�Ծ�����֤��ͻ��ԣ�����������ʱȶ�ʱҲӦ���մ�ԭ���Ӧ�ظ������顣
	int mark_saving = 0;
	for (int i = 0; i < list_result.size();i++) {
		spectra& pep = list_result[i];//����򵥵ı��������̣���߿ɶ���
		//if (pep.marker != mark_saving) {
			if (pep.is_mut == false) {
				if (pimap.find(pep.prot) != pimap.end()) {//�ҵ��˵Ļ�
					vector<int> j = pimap[pep.prot];
					//pos_saving = j;//���ݸ�pos_saving
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
						//�Ҳ����Ļ���������Ȼ�����ѭ��
					}
			}
			else {//��������
				if (pimap.find(pep.prot)!=pimap.end()) {//���Һ�ģ���Ч
					vector<int> j = pimap[pep.prot];
					//������openresearch����ʱ�����ð������
					///mut_pep_inform temp2 = pepmutation(pep,intri,table);//ͻ������
					for (int i = 0; i < j.size(); i++) {//ѭ��ͻ�䵰�������в��ȶ�
						string temp = list_pro[j[i]].hseq;
						temp = temp.replace(list_pro[j[i]].pos, 1, 1, list_pro[j[i]].mutataa);
						///auto pos_find = temp.find(temp2.mutpep);
					//����������ͻ������ʱ�����ð������
						auto pos_find = temp.find(pep.seq);
						///bool access = false;//�˱��������ж��Ƿ�������һ���Ķ�ͻ��λ���뵰����ͻ��λ���غ�
						///for (int i = 0; i < temp2.size; i++) {//�˴�pos_mut��result���Ķ�ͻ��λ�������(0��ʼ)����ƥ����ʼλ�����꣨pos_find����ͬ����0��ʼ�����Ӧ���ڵ�����������ͻ��λ�����꣨0��ʼ��
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
	/////���
	string outname = "openresearchɸѡ��.txt";
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

