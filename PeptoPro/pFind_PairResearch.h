#pragma once
#include<string>
#include<vector>
using namespace std;
struct spectra
{
	string file_name;
	int scan_no;
	double exp_mh;
	int charge;
	double q_value;
	string seq;
	double calc_mh;
	double mass_shift;
	double raw_score;
	string final_score;
	string modi;
	int spec;
	string prot;
	string posi;
	string label;
	string targe;
	int mc_sites;
	double afm_shift;
	int others;
	int marker = 0;
	bool is_modi = false;
};

ostream& operator<<(ostream& os, const spectra& s);
istream& operator>>(istream& os, vector<spectra>& list_result);
void mark(vector<spectra>& list);
bool sortbyleg(spectra& a, spectra& b);
bool sortbymarker(spectra& a, spectra& b);