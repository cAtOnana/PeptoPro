#include<vector>
#include"reflect.h"
#include"pFind_PairResearch.h"
using namespace std;

int main(int argc, char* argv[]) {//result mrna��ͬ��ͻ�� ���ձ� ����������fasta
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
}


ostream & operator<<(ostream & os, spectra & s)
{
	os << s.file_name << "	" << s.scan_no << "	" << s.exp_mh << "	" << s.charge << "	" << s.q_value << "	" << s.seq << "	"
		<< s.calc_mh << "	" << s.mass_shift << "	" << s.raw_score << "	" << s.final_score << "	" << s.modi << "	" << s.spec
		<< "	" << s.prot << "	" << s.posi << "	" << s.label << "	" << s.targe << "	" << s.mc_sites << "	" << s.afm_shift << "	"
		<< s.others;
	return os;
}