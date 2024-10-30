/*

This code was prepared to generate the cosine similarity scores and min-max
indices using cosine similarity scores presented in the manuscript 
"The min-max test: an objective method for discriminating mass spectra" by
Moorthy and Sisco (2021).

This code was prepared for research purposes and is licenced as NIST software
(see license file in parent folder).

The code was prepared by Arun Moorthy (arun.moorthy@nist.gov).
Questions or comments can be directed towards him.

The code was developed and tested in unix environments (MacOSX and Ubuntu 18.04)
Date: 2021/05/28
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <numeric>

using namespace std;

const int MAXRELINTENSITY = 1;
const int MAXNUMPEAKS = 1000;  // should only be considering small molecules less that 1000 Daltons
const int DEBUG = 0;


// Spectrum Class Definition and function definitions
class Spectrum
{
	string spectrumFile;
	int numpeaks;
	vector<double> mz;
	vector<double> ab;
	public:
	Spectrum(string);
	~Spectrum();
	string getFN() {return spectrumFile;}
	vector<double> getMZ() {return mz;}
	vector<double> getAB() {return ab;}
	int getNP() {return numpeaks;}
	vector<double> lowResVec(double exp1, double exp2);

};

Spectrum::Spectrum(string spectrumFile)
{
	this -> spectrumFile = spectrumFile;

	int 	i = 0, j = 0, k = 0; // indices for future use
	int 	numpeaks=0;

	string line;
	char letter;
	char commentChar = '#';
	string s_number;
	double d_number;

	vector<double> mz;
	vector<double> ab;

	ifstream spectrum(spectrumFile);
	if (!spectrum){
		cout << "Error reading query data file.\n";
		// return -1;
	} else {
		for(i=0;!spectrum.eof();i++)
		{
			getline(spectrum,line);
			letter = line[0];
			if(letter != commentChar){
				numpeaks++;
				k = 0;
				for(j=0;j<line.length();j++){
					if(line[j] == ','){

						if(k==1){
							d_number = stod(s_number);
							// cout << "mz: " << d_number << "\t";
							mz.push_back(d_number);
							s_number.clear();
						}

						k++;
						s_number.clear();

					} else if (j==(line.length()-1)) {
						d_number = stod(s_number);
						// cout << "ab: " << d_number << "\n";
						ab.push_back(d_number);
						s_number.clear();
					} else {
						s_number += line[j];
					}

				}
				s_number.clear();
			}

		}
	}

	spectrum.close();


	if(DEBUG==1){
		cout << "Iterations: " << i << endl;
		cout << "Num Peaks: " << numpeaks << endl << endl;

		for(i=0;i<mz.size();i++){
			cout << mz[i] << "\t" << ab[i] << endl;
		}
	}

	this -> numpeaks = numpeaks;
	this -> mz = mz;
	this -> ab = ab;

	// return 0;
}

//Spectrum::~Spectrum() { cout << spectrumFile <<" spectrum destroyed.\n";}
Spectrum::~Spectrum() { }

vector<double> Spectrum::lowResVec(double exp1, double exp2)
{

	int i = 0;
	int int_mz=0;

	vector<double> weightedAB(ab.size(),0.0);
	for(i = 0;i<ab.size();i++){
		// cout << round(mz[i]) << "\t" << ab[i] << "\t";
		weightedAB[i] = pow(round(mz[i]),exp1) * pow(ab[i],exp2);
		// cout << weightedAB[i] << endl;
	}


	double maxab = *max_element(weightedAB.begin(),weightedAB.end());
	vector<double> result(MAXNUMPEAKS,0.0);

	for(i=0;i<mz.size();i++){
		int_mz = (int) round(mz[i]);
		result[int_mz] += MAXRELINTENSITY*weightedAB[i]/maxab;
	}

	if(DEBUG==1){
		for(i = 0; i<result.size();i++){
			if(result[i]!=0){
				cout << i << "\t" << result[i] << endl;
			}
		}
	}

	return result;
}


// Global functions

vector<string> asm_randomsample(vector<string> x, int n)
{
	vector<string> result;
	int irand;
	string strand;

	int i=0;
	int it=0;

			srand (time(NULL));

	while(i<n){

		it++;
		irand = rand() % 10 ;
		strand = x[irand];
			if (find(result.begin(), result.end(), strand) == result.end()){
					result.push_back(strand);
					i++;
			}

	}

	return result;
}

double cosine_similarity(const vector<double> &x, const vector<double> &y, const int &specLength)
{
	// function as described in Moorthy & Kearsley 2021.
	double num = 0.0, den1 = 0.0, den2 = 0.0, result = 0.0;
	int i = 0;

	for(i =0; i<specLength; i++){
		// cout << x[i] << "\t" << y[i] << endl;
		num  += x[i] * y[i];
        den1 += x[i] * x[i];
		den2 += y[i] * y[i];
	}

	result = num / (sqrt(den1) * sqrt(den2));
	return result;
}

double simple_MF(const vector<double> &x, const vector<double> &y, const int &specLength)
{
	// function as described in Moorthy & Kearsley 2021.

	double 	result = 0.0;

	result = cosine_similarity(x,y,specLength);
	result = MAXRELINTENSITY*result*result;
	return result;
}

double identity_MF(const vector<double> &x, const vector<double> &y, const int &specLength)
{
	// function as described in Moorthy & Kearsley 2021.

	int i = 0;
	vector<double> r;
	double rnum=0.0;
	double test = 0.0;

	// compute peak pair ratio
	for(i=1;i<specLength;i++){


		test = x[i]*x[i-1]*y[i]*y[i-1]; // if any of these elements are zero, advance loop

		if(test>1e-8){
			rnum = x[i]*y[i-1]/(x[i-1]*y[i]);
			r.push_back(rnum);
		} else {
			r.push_back(0);
		}
	}

	double Fnum=0.0;
	double Fden=0.0;
	double F = 0.0;
	double m1=0;

	for(i=0;i<r.size();i++){
		if(r[i]>1e-8){
			m1++;
			Fnum += i * min(r[i],1/r[i]);
			Fden += i;
		}
	}

	F = Fnum/(Fden+1e-8);  // prevents nan when there are no adjacent peaks

	double m2=0;
	for(i=0;i<specLength;i++){
		test = x[i]*y[i];
		if(test>1e-8){
			m2++;
		}
	}

	// put peak pair ratio information together with simple similarity
	double result = 0.0;
	double epsilon2 = simple_MF(x,y,specLength);

	result = MAXRELINTENSITY*(m1*F + m2*epsilon2/MAXRELINTENSITY)/(m1+m2+1e-8);
	// result = (m1*F + m2*epsilon2/MAXRELINTENSITY)/(m1+m2);
	// result = (m1*F + m2*epsilon2)/(m1+m2);
	//	cout << "\tm1: " << m1 << endl;
	//	cout << "\tm2: " << m2 << endl;
	//	cout << "\tF: " << F << endl;
	//	cout << "\tepsilon2/c: " << epsilon2/MAXRELINTENSITY << endl;

	return result;

}

double MaxNorm(const vector<double> &x, const vector<double> &y, const int &specLength)
{
	// function as described in Moorthy & Kearsley 2021.
	double result = 0.0;
	double value = 0.0;
	int i = 0;

	for(i=0;i<specLength;i++){
		value = x[i] - y[i];
		value = fabs(value);
		result = max(result,value);
	}

	return result;

}

double L1Norm(const vector<double> &x, const vector<double> &y, const int &specLength)
{
	// function as described in Moorthy & Kearsley 2021.
	double result = 0.0;
	double value = 0.0;
	int i = 0;

	for(i=0;i<specLength;i++){
		value = x[i] - y[i];
		value = fabs(value);
		result += value;
	}

	return result;

}

double L2Norm(const vector<double> &x, const vector<double> &y, const int &specLength)
{
	// function as described in Moorthy & Kearsley 2021.
	double result = 0.0;
	double value = 0.0;
	int i = 0;

	for(i=0;i<specLength;i++){
		value = (x[i] - y[i]) * (x[i]-y[i]) ;
		result += value;
	}

	result = sqrt(result);
	return result;

}

double asm_minmaxTest(double minS11, double minS22, double maxS12)
{
	double result = 0.0;

	 if(minS11<minS22){
		 result = minS11 - maxS12;
	 } else {
		 result = minS22 - maxS12;
	 }

	 result = 1 - max(result,0.0); // transforms to scale between 0 and 1

	 //result = (1.0 - result)/2.0; // transform to scale between 0 and 1.
	 //cout << minS11 << " " << minS22 << " " << maxS12 << " " << result << endl;

	return result;
}

vector<double> asm_statistics(vector<double> x)
{
	vector<double> result;
	int i=0;
	int sizeX = x.size();

	sort(x.begin(),x.end());


	double minx = *min_element(x.begin(),x.end());
	//cout << minx << endl;
	result.push_back(minx);


	double medianx=0.0;
	 if (sizeX % 2 != 0)
        medianx = x[sizeX / 2];
 	else{
 		medianx = (x[(sizeX - 1) / 2] + x[sizeX / 2]) / 2.0;
 	}
 	//cout << medianx << endl;
 	result.push_back(medianx);



	double maxx = *max_element(x.begin(),x.end());
	//cout << maxx << endl;
	result.push_back(maxx);

	double meanx=0.0;

	for(i=0;i<sizeX;i++){
		meanx += x[i];
	}
	meanx = meanx/(double) sizeX;
	//cout << meanx << endl;
	result.push_back(meanx);

	double varx=0.0;
	double stdevx=0.0;

	for(i=0;i<sizeX;i++){
		varx += pow((x[i]-meanx),2);
	}
	varx = varx/(double) (sizeX - 1);
	stdevx = sqrt(varx);
	//cout << stdevx << endl << endl;
	result.push_back(stdevx);

	return result;

}

vector<double> asm_ROC(vector<double> &sortedX, vector<int> &sortedClass)
{
	vector<double> result;

	int i;
	int index;

	vector<double> thresholds;
	double increment=0;
	for(i=0;i<1000;i++){
		increment += 0.001;
		if((sortedX[0]-increment) < 1e-8){
			if((sortedX.back()-increment) > -1e-8){
				thresholds.push_back(increment);
			}
		}
	}

	double tn=0;
	double fn=0;
	double tp=0;
	double fp=0;

	double accuracy;
	double accuracy_opt=0;
	double recall;
	double specificity;
	double Fpr;
	double precision;
	double Fmeasure;

	//int i = 0;
	int i_opt = 0;

	for(i=0;i<thresholds.size();i++){

		index = upper_bound(sortedX.begin(),sortedX.end(), thresholds[i]) - sortedX.begin()-1;

		tn = count(sortedClass.begin(),sortedClass.begin()+index,0);
		fn = count(sortedClass.begin(),sortedClass.begin()+index,1);
		fp = count(sortedClass.begin()+index,sortedClass.end(),0);
		tp = count(sortedClass.begin()+index,sortedClass.end(),1);

     /*
		cout << "Threshold: " << sortedX[i];
		cout << " TP: " << tp ;
		cout << " FP: " << fp;
		cout << " TN: " << tn;
		cout << " FN: " << fn << endl;
     */

		accuracy = (tp+tn)/(tp+tn+fn+fp);
		recall = tp/(tp+fn);
		specificity = tn/(tn+fp);
		Fpr = fp/(tn+fp);
		precision = tp/(tp+fp);
		Fmeasure = 2.0/(1/precision + 1/recall);

		if(accuracy > accuracy_opt){
			i_opt = index;
			accuracy_opt = accuracy;
		}

		tn=0;fn=0;tp=0;fp=0;index=0;
	}

	i = i_opt;
	tn = count(sortedClass.begin(),sortedClass.begin()+i,0);
	fn = count(sortedClass.begin(),sortedClass.begin()+i,1);
	fp = count(sortedClass.begin()+i,sortedClass.end(),0);
	tp = count(sortedClass.begin()+i,sortedClass.end(),1);

	accuracy = (tp+tn)/(tp+tn+fn+fp);
	recall = tp/(tp+fn);
	specificity = tn/(tn+fp);
	Fpr = fp/(tn+fp);
	precision = tp/(tp+fp);
	Fmeasure = 2.0/(1/precision + 1/recall);

	//cout << endl;
	//cout << "Threshold (opt accuracy): " << sortedX[i_opt] << endl;
	result.push_back(sortedX[i_opt]);
	//cout << "Accuracy: " << accuracy << endl;
	result.push_back(accuracy);
	//cout << "Recall: " << recall << endl;
	result.push_back(recall);
	//cout << "Specificity: " << specificity << endl;
	result.push_back(specificity);
	//cout << "FPR: " << Fpr << endl;
	result.push_back(Fpr);
	//cout << "Precision: " << precision << endl;
	result.push_back(precision);
	//cout << "F-measure: " << Fmeasure << endl;
	result.push_back(Fmeasure);

	//cout << endl << endl;

	return result;
}

void MF_report(const vector<double> &MF, const vector<int> &IsRep, const double &duration, ofstream &writer)
{
	auto t1 = chrono::high_resolution_clock::now();

  writer << "Summary of Match Factor threshold study\n";
	int i = 0;

  time_t now = time(0);
  char* dt = ctime(&now);
  writer << dt << endl << endl;

	// standard statistics
	vector<double> TruePosMF;
	vector<double> TrueNegMF;

	for(i=0;i<MF.size();i++){
		if(IsRep[i]==1){
			TruePosMF.push_back(MF[i]);
		} else {
			TrueNegMF.push_back(MF[i]);
		}
	}


	vector<double> TruePosStats=asm_statistics(TruePosMF);
	vector<double> TrueNegStats=asm_statistics(TrueNegMF);

	writer << "Frequency: " << TruePosMF.size() << " (rep); " << TrueNegMF.size() << " (other)" << endl;
	writer << "Min MF: "	<< TruePosStats[0] << " (rep); " << TrueNegStats[0] << " (other)" << endl;
	writer << "Median MF: "	<< TruePosStats[1] << " (rep); " << TrueNegStats[1] << " (other)" << endl;
	writer << "Max MF: "	<< TruePosStats[2] << " (rep); " << TrueNegStats[2] << " (other)" << endl;
	writer << "Mean MF: "	<< TruePosStats[3] << " (rep); " << TrueNegStats[3] << " (other)" << endl;
	writer << "SD MF: "	<< TruePosStats[4] << " (rep); " << TrueNegStats[4] << " (other)" << endl;
  writer << endl << endl;
	// ROC statistics
	vector<int> indices(IsRep.size());
	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(),[&](int A, int B) -> bool {return MF[A] < MF[B];});

	vector<double> sortedMF;
	vector<int> sortedIsRep;

	for(i=0;i<indices.size();i++){
		sortedMF.push_back(MF[indices[i]]);
		sortedIsRep.push_back(IsRep[indices[i]]);
	}

  vector<double> ROC_results = asm_ROC(sortedMF,sortedIsRep);

  writer << "Threshold (opt accuracy): " << ROC_results[0] << endl;
	writer << "Accuracy: " << ROC_results[1] << endl;
  writer << "Recall: " << ROC_results[2] << endl;
	writer << "Specificity: " << ROC_results[3] << endl;
	writer << "FPR: " << ROC_results[4] << endl;
	writer << "Precision: " << ROC_results[5] << endl;
	writer << "F-measure: " << ROC_results[6] << endl;

  writer << endl << endl;

	auto t2 = chrono::high_resolution_clock::now();

  auto duration2 = chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

	writer << "MF compute time: " << duration/1000000 << " seconds.\n";
	writer << "ROC and Analysis time: " << duration2/1000000 << " seconds.\n";
}

void MinMax_report(const vector<double> &MF, const vector<int> &IsRep, const double &duration, ofstream &writer)
{
	auto t1 = chrono::high_resolution_clock::now();

  writer << "Summary of MinMax threshold study\n";
	int i = 0;

  time_t now = time(0);
  char* dt = ctime(&now);
  writer << dt << endl << endl;

	// standard statistics
	vector<double> TruePosMF;
	vector<double> TrueNegMF;

	for(i=0;i<MF.size();i++){
		if(IsRep[i]==1){
			TruePosMF.push_back(MF[i]);
		} else {
			TrueNegMF.push_back(MF[i]);
		}
	}


	vector<double> TruePosStats=asm_statistics(TruePosMF);
	vector<double> TrueNegStats=asm_statistics(TrueNegMF);

	writer << "Frequency: " << TruePosMF.size() << " (rep); " << TrueNegMF.size() << " (other)" << endl;
	writer << "Min Delta12: "	<< TruePosStats[0] << " (rep); " << TrueNegStats[0] << " (other)" << endl;
	writer << "Median Delta12: "	<< TruePosStats[1] << " (rep); " << TrueNegStats[1] << " (other)" << endl;
	writer << "Max Delta12: "	<< TruePosStats[2] << " (rep); " << TrueNegStats[2] << " (other)" << endl;
	writer << "Mean Delta12: "	<< TruePosStats[3] << " (rep); " << TrueNegStats[3] << " (other)" << endl;
	writer << "SD Delta12: "	<< TruePosStats[4] << " (rep); " << TrueNegStats[4] << " (other)" << endl;
  writer << endl << endl;

	// ROC statistics
	vector<int> indices(IsRep.size());
	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(),[&](int A, int B) -> bool {return MF[A] < MF[B];});

	vector<double> sortedMF;
	vector<int> sortedIsRep;

	for(i=0;i<indices.size();i++){
		sortedMF.push_back(MF[indices[i]]);
		sortedIsRep.push_back(IsRep[indices[i]]);
	}

  vector<double> ROC_results = asm_ROC(sortedMF,sortedIsRep);
  writer << "Threshold (opt accuracy): " << ROC_results[0] << endl;
	writer << "Accuracy: " << ROC_results[1] << endl;
  writer << "Recall: " << ROC_results[2] << endl;
	writer << "Specificity: " << ROC_results[3] << endl;
	writer << "FPR: " << ROC_results[4] << endl;
	writer << "Precision: " << ROC_results[5] << endl;
	writer << "F-measure: " << ROC_results[6] << endl;

  writer << endl << endl;

	auto t2 = chrono::high_resolution_clock::now();

  auto duration2 = chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

	writer << "MinMax compute time: " << duration/1000000 << " seconds.\n";
	writer << "ROC and Analysis time: " << duration2/1000000 << " seconds.\n";
}

void BuildListOfFiles(string List, vector<string> &c1)
{
	int i, j, k, ii, jj, kk; // for index
	string line;
	char letter;
	char commentChar = '#';
	string s_number;

	ifstream compareList(List);
	if (!compareList){
		cout << "Error reading list of files.\n";
	} else {
		for(i=0;!compareList.eof();i++)
		{
			getline(compareList,line);
			letter = line[0];
			if(letter != commentChar){
				k = 0;
				s_number.clear();
				for(j=0;j<line.length();j++){
					if(line[j] == ','){

						if(k==0){
							s_number += "/";
						}

						if(k==1){
							//cout << s_number << endl;
							c1.push_back(s_number);
							s_number.clear();
						}

						k++;

					} else {
						s_number += line[j];
					}

				}

			}

		}

	}

	compareList.close();





}

int main()
{
	int flag = 0;
	int i, j, k, ii, jj, kk; // for indexing
	double mf;

	// /*
	vector<string> ListOfFiles;
	BuildListOfFiles("CodeNames_withFormula.csv",ListOfFiles);
	int numFiles = ListOfFiles.size();
	vector<string> replicate {"01","02","03","04","05","06","07","08","09","10"};
	int numReplicates = replicate.size();
	string queryname;
	string refname;
	double totIter = numFiles*numReplicates*numFiles*numReplicates;
	double iter = 0;

	// Match Factor statistic (standard analysis trying to threshold by match factor)
	vector<double> MF;
	vector<int> IsRep;

	ofstream mfoutput("MF_fullResults.txt");

	auto t1 = chrono::high_resolution_clock::now();
	for(i=0;i<numFiles;i++){
		for(ii=0;ii<numReplicates;ii++){

			queryname = ListOfFiles[i] + "-" + replicate[ii] + ".csv";
			Spectrum query(queryname);
			vector<double> x = query.lowResVec(0,0.5);

			for(j=0;j<numFiles;j++){
				for(jj=0;jj<numReplicates;jj++){

					iter++;
					cout << iter/totIter << "\t" << queryname << "\t" << refname << "\n";

					if((i!=j) || (ii!=jj)){

						refname = ListOfFiles[j] + "-" + replicate[jj] + ".csv";

						mfoutput << ListOfFiles[i] << "\t" << ListOfFiles[j] << "\t";

						Spectrum library(refname);
						vector<double> y = library.lowResVec(0,0.5);

									mf = cosine_similarity(x,y,MAXNUMPEAKS);

												//if(fabs(1.0-mf)>1e-8){
																 MF.push_back(mf);
																 mfoutput << mf << "\t";
																				if(i==j){
																					IsRep.push_back(1);
																					mfoutput << "1" << endl;
																				} else {
																					IsRep.push_back(0);
																					mfoutput << "0" << endl;
																				}
												//}
					}
				}
			}
		}
	}
	mfoutput.close();

	auto t2 = chrono::high_resolution_clock::now();

	auto duration = chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

	ofstream writer("MFReport.txt");
	MF_report(MF,IsRep,duration,writer);
	writer.close();


  // DELTA12 statistic (min-max test)
	vector<double> DELTA12;
	IsRep.clear();
	double delta12;
	double minS11=MAXRELINTENSITY;
	double minS22=MAXRELINTENSITY;
	double maxS12=0.0;
	int numRepsInTest = 3; //replicate.size()/2; // maximum value
	int numTests = 100; // increase to 100 to compare with MF threshold study

	vector<string> repsToStudy;
	vector<string> repsToStudyQ;
	vector<string> repsToStudyL;

	ofstream minmaxoutput("MinMax_fullResults.txt");

	totIter = numFiles*numFiles*numTests;
	iter = 0;
	t1 = chrono::high_resolution_clock::now();
	for(ii = 0;ii < numTests; ii++){
	for(i=0;i<numFiles;i++){
			for(j=0;j<numFiles;j++){

						iter++;
						cout << iter/totIter << "\n";

						repsToStudy = asm_randomsample(replicate,numRepsInTest*2);
						repsToStudyQ.insert(repsToStudyQ.begin(),repsToStudy.begin(),repsToStudy.begin()+numRepsInTest);
						repsToStudyL.insert(repsToStudyL.begin(),repsToStudy.begin()+numRepsInTest,repsToStudy.end());
						//cout << repsToStudy.size() << " " << repsToStudyQ.size() << " " << repsToStudyL.size() << endl;

						// compute minS11
						mf = 1.0;
						for(k=0;k < repsToStudyQ.size();k++){
							queryname = ListOfFiles[i] + "-" + repsToStudyQ[k] + ".csv";
							Spectrum query(queryname);
							vector<double> x = query.lowResVec(0,0.5);
							for(kk = 0;kk<repsToStudyQ.size();kk++){
								if(kk != k){
									refname = ListOfFiles[i] + "-" + repsToStudyQ[kk] + ".csv";
									Spectrum library(refname);
									vector<double> y = library.lowResVec(0,0.5);

												mf = min(mf,cosine_similarity(x,y,MAXNUMPEAKS));

								}
							}
						}
						minS11 = mf;

						// compute minS22
						mf = 1.0;
						for(k=0;k < repsToStudyL.size();k++){
							queryname = ListOfFiles[j]+"-"+repsToStudyL[k]+".csv";
							Spectrum query(queryname);
							vector<double> x = query.lowResVec(0,0.5);
							for(kk = 0;kk<repsToStudyL.size();kk++){
								if(kk != k){
									refname = ListOfFiles[j]+"-"+repsToStudyL[kk]+".csv";
									Spectrum library(refname);
									vector<double> y = library.lowResVec(0,0.5);

												mf = min(mf,cosine_similarity(x,y,MAXNUMPEAKS));

								}
							}
						}
						minS22 = mf;

						// compute maxS12
						mf = 0.0;
						for(k=0;k < repsToStudyQ.size();k++){
							queryname = ListOfFiles[i]+"-"+repsToStudyQ[k]+".csv";
							Spectrum query(queryname);
							vector<double> x = query.lowResVec(0,0.5);
							for(kk = 0;kk<repsToStudyL.size();kk++){


									refname = ListOfFiles[j]+"-"+repsToStudyL[kk]+".csv";
									Spectrum library(refname);
									vector<double> y = library.lowResVec(0,0.5);

												mf = max(mf,cosine_similarity(x,y,MAXNUMPEAKS));
												mf = max(mf,cosine_similarity(y,x,MAXNUMPEAKS));
												// cout << k << " " << kk << " " << mf << endl;

							}
						}
						maxS12 = mf;

						minmaxoutput << ListOfFiles[i] << "\t" << ListOfFiles[j] << "\t";
						minmaxoutput << minS11 << "\t" << minS22 << "\t" << maxS12 << "\t";
						delta12 = asm_minmaxTest(minS11,minS22,maxS12);
						minmaxoutput << delta12 << "\t";

						DELTA12.push_back(delta12);
						if(i==j){
							IsRep.push_back(1);
							minmaxoutput<<"1"<<endl;
						} else {
							IsRep.push_back(0);
							minmaxoutput<<"0"<<endl;
						}

						minS11=MAXRELINTENSITY;
						minS22=MAXRELINTENSITY;
						maxS12=0.0;
						mf = MAXRELINTENSITY;
						repsToStudy.clear();
						repsToStudyQ.clear();
						repsToStudyL.clear();

			}

  }
	}
	cout << iter << endl;
	t2 = chrono::high_resolution_clock::now();

  duration = chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();



  ofstream writer2("MinMaxReport.txt");
  MinMax_report(DELTA12,IsRep,duration,writer2);
  writer2.close();



return 0;
}
