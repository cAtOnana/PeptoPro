#include"reflect.h"
#include"pFind_PairResearch.h"
#include<vector>
#include<algorithm>
//result�ļ����Ķ�����nm�ű�ǵģ���ͬ��ͻ���ļ��еĵ�������������ensp�ű�ǵģ���ͨ������հ����Ķ��ձ���ӳ�佫���߶�Ӧ��������бȶԣ���ͻ����Ķ�ֱ�ӱȶԼ��ɣ���ͻ����Ķ����ڶ�Ӧ�����������н���ͻ��ת��Ȼ��ȶ�
using namespace std;
int argc = 5;
char* argv[5] = { "aaa","result.txt","��ͬ��ͻ��.txt","ENSP-NM.txt","Homo_sapiens.GRCh38.pep.all_cnew.fa" };
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
	int saving = 0;
	///��¼ƥ�䵽��list_pro�����꣬������һ��ƥ��
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
	for (int i = 0; i < list_result.size();) {//���·��ڷ�֧��
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
