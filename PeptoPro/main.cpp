#include"reflect.h"
#include"pFind_PairResearch.h"
#include<vector>
#include<unordered_map>
//result�ļ����Ķ�����nm�ű�ǵģ���ͬ��ͻ���ļ��еĵ�������������ensp�ű�ǵģ���ͨ������հ����Ķ��ձ���ӳ�佫���߶�Ӧ��������бȶԣ���ͻ����Ķ�ֱ�ӱȶԼ��ɣ���ͻ����Ķ����ڶ�Ӧ�����������н���ͻ��ת��Ȼ��ȶ�
using namespace std;
int argc = 5;
char* argv[5] = { "aaa","result.txt","��ͬ��ͻ��.txt","ENSP-NM.txt","Homo_sapiens.GRCh38.pep.all_cnew.fa" };
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
		return a1.find(a2)!=string::npos;//���Դ˴���������ж϶��ǰ����ж�ʱ�Ƿ����
	}
};
typedef unordered_map<string, int, pro_hash, pro_compare> ProIndexMap;
int main(){//int argc, char* argv[]) {//result mrna��ͬ��ͻ�� ���ձ� ����������fasta
	char* logname = "PeptoPro.log";
	ofstream log(logname);
	if (argc != 5) {
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
	vector<spectra> list_result;
	inres >> list_result;
	int markcount =mark(list_result);
	for (int i = 0; i < list_result.size();i++)
		log << list_result[i];

	vector<pro> list_pro;
	inmrn >> list_pro;
	fillhseq(inref, list_pro);
	//������ձ�
	fillnm(incom, list_pro);
	//��ʼ�ȶ�
	bool* cout_modi = new bool[markcount] {false};
	bool* cout_no = new bool[markcount] {false};
	///�����������ڼ�¼ĳһmark�Ŷ�Ӧ���������ܷ������ֻ�����߽�Ϊtrue����ĳһ����������һ�Ծ�����֤��ͻ��ԣ�����������ʱȶ�ʱҲӦ���մ�ԭ���Ӧ�ظ������顣
	int pos_saving = 0, mark_saving = 0;
	///pos_saving��¼ƥ�䵽��list_pro�����꣬������һ��ƥ�䣻mark_saving���ڼ�¼��ǰ���
	///mark_saving��¼��ǰ����ִ��ƥ����������������飻ͬʱ�����ж�pos_saving��ʹ����񣬽���߽����⣨��Ϊpos_saving�ǽ���ֱ�ӱȶ����̵ģ�
	ProIndexMap pimap;
	for (int i = 0; i < list_pro.size(); i++) {
		pimap[list_pro[i].nm] = i;
	}
	for (int i = 0; i < list_result.size();) {//���·��ڷ�֧ĩ
		spectra& pep = list_result[i];//����򵥵ı��������̣���߿ɶ���
		if (pep.marker != mark_saving) {
			if (pep.is_modi == false) {
				if (pimap.find(pep.prot) != pimap.end()) {//�ҵ��˵Ļ�
					int j = pimap[pep.prot];
					if (list_pro[j].nm.find(pep.prot) != string::npos) {
						pos_saving = j;//���ݸ�pos_saving
						mark_saving =pep.marker;
						cout_no[mark_saving] = true;
						i++;//update
						break;
					}
				}
					else {//�Ҳ����Ļ���������Ȼ�����ѭ��
						vector<spectra>::iterator itor = list_result.begin() + i;
						list_result.erase(itor);//erase & update
					}
			}
			else {//��������
				if (pimap.find(pep.prot)!=pimap.end()) {//���Һ�ģ���Ч
					int j = pimap[pep.prot];
					string temp = list_pro[j].hseq;
					temp = temp.replace(list_pro[j].pos, 1, 1, list_pro[j].mutataa);
					string temp2 =pep.prot;//A ?��Ŀ���������ж��ͻ�䣬��ô�죿   �𣺺ϲ�ȫ��ͻ��
					auto pos_find = temp.find(temp2);
					if (pos_find != string::npos && pos_find + pos_mut == list_pro[j].pos) {//�˴�pos_mut��result���Ķ�ͻ��λ������꣬��ƥ����ʼλ�����꣨pos_find�����Ӧ���ڵ�����������ͻ��λ������
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
		else {
			if (pep.is_modi == false) {
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
				string temp = list_pro[pos_saving].hseq;
				temp = temp.replace(list_pro[pos_saving].pos, 1, 1, list_pro[pos_saving].mutataa);//����A
				string temp2 = pep.prot;
				//�˴���temp2��Ϊͻ������
				if (list_pro[pos_saving].nm.find(pep.prot) != string::npos) {
					auto pos_find = temp.find(temp2);
					if (pos_find != string::npos && pos_find + pos_mut == list_pro[j].pos) {//�˴�pos_mut��result���Ķ�ͻ��λ������꣬��ƥ����ʼλ�����꣨pos_find�����Ӧ���ڵ�����������ͻ��λ������
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

