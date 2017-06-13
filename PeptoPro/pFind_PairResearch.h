#pragma once
#include<string>
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
};

ostream& operator<<(ostream& os, spectra& s);
istream& operator>>(istream& os, spectra& s);