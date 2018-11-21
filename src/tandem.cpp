/*
 Copyright (C) 2003-2013 Ronald C Beavis, all rights reserved
 X! tandem 
 This software is a component of the X! proteomics software
 development project

Use of this software governed by the Artistic license, as reproduced here:

The Artistic License for all X! software, binaries and documentation

Preamble
The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some 
semblance of artistic control over the development of the package, 
while giving the users of the package the right to use and distribute 
the Package in a more-or-less customary fashion, plus the right to 
make reasonable modifications. 

Definitions
"Package" refers to the collection of files distributed by the Copyright 
	Holder, and derivatives of that collection of files created through 
	textual modification. 

"Standard Version" refers to such a Package if it has not been modified, 
	or has been modified in accordance with the wishes of the Copyright 
	Holder as specified below. 

"Copyright Holder" is whoever is named in the copyright or copyrights 
	for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of 
	media cost, duplication charges, time of people involved, and so on. 
	(You will not be required to justify it to the Copyright Holder, but 
	only to the computing community at large as a market that must bear 
	the fee.) 

"Freely Available" means that no fee is charged for the item itself, 
	though there may be fees involved in handling the item. It also means 
	that recipients of the item may redistribute it under the same
	conditions they received it. 

1. You may make and give away verbatim copies of the source form of the 
Standard Version of this Package without restriction, provided that 
you duplicate all of the original copyright notices and associated 
disclaimers. 

2. You may apply bug fixes, portability fixes and other modifications 
derived from the Public Domain or from the Copyright Holder. A 
Package modified in such a way shall still be considered the Standard 
Version. 

3. You may otherwise modify your copy of this Package in any way, provided 
that you insert a prominent notice in each changed file stating how and 
when you changed that file, and provided that you do at least ONE of the 
following: 

a.	place your modifications in the Public Domain or otherwise make them 
	Freely Available, such as by posting said modifications to Usenet 
	or an equivalent medium, or placing the modifications on a major 
	archive site such as uunet.uu.net, or by allowing the Copyright Holder 
	to include your modifications in the Standard Version of the Package. 
b.	use the modified Package only within your corporation or organization. 
c.	rename any non-standard executables so the names do not conflict 
	with standard executables, which must also be provided, and provide 
	a separate manual page for each non-standard executable that clearly 
	documents how it differs from the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

4. You may distribute the programs of this Package in object code or 
executable form, provided that you do at least ONE of the following: 

a.	distribute a Standard Version of the executables and library files, 
	together with instructions (in the manual page or equivalent) on 
	where to get the Standard Version. 
b.	accompany the distribution with the machine-readable source of the 
	Package with your modifications. 
c.	give non-standard executables non-standard names, and clearly 
	document the differences in manual pages (or equivalent), together 
	with instructions on where to get the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

5. You may charge a reasonable copying fee for any distribution of 
this Package. You may charge any fee you choose for support of 
this Package. You may not charge a fee for this Package itself. 
However, you may distribute this Package in aggregate with other 
(possibly commercial) programs as part of a larger (possibly 
commercial) software distribution provided that you do not a
dvertise this Package as a product of your own. You may embed this 
Package's interpreter within an executable of yours (by linking); 
this shall be construed as a mere form of aggregation, provided that 
the complete Standard Version of the interpreter is so embedded. 

6. The scripts and library files supplied as input to or produced as 
output from the programs of this Package do not automatically fall 
under the copyright of this Package, but belong to whomever generated 
them, and may be sold commercially, and may be aggregated with this 
Package. If such scripts or library files are aggregated with this 
Package via the so-called "undump" or "unexec" methods of producing 
a binary executable image, then distribution of such an image shall 
neither be construed as a distribution of this Package nor shall it 
fall under the restrictions of Paragraphs 3 and 4, provided that you 
do not represent such an executable image as a Standard Version of 
this Package. 

7. C subroutines (or comparably compiled subroutines in other languages) 
supplied by you and linked into this Package in order to emulate 
subroutines and variables of the language defined by this Package 
shall not be considered part of this Package, but are the equivalent 
of input as in Paragraph 6, provided these subroutines do not change 
the language in any way that would cause it to fail the regression 
tests for the language. 

8. Aggregation of this Package with a commercial distribution is always 
permitted provided that the use of this Package is embedded; that is, 
when no overt attempt is made to make this Package's interfaces visible 
to the end user of the commercial distribution. Such use shall not be 
construed as a distribution of this Package. 

9. The name of the Copyright Holder may not be used to endorse or promote 
products derived from this software without specific prior written permission. 

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE. 

The End 
*/

// File version: 2013-04-01

/*
	tandem.cpp provides a command line interface to the main processing class, mprocess.
	A single command line parameter is accepted, which is the file name for an XML class
	containing the parameters for performing the protein modelling session. The input file 
	also contains the name of an output file, which will contain the output from the session.
*/

#include "stdafx.h"
#include <sys/timeb.h>
#include <sys/time.h>
#include <ctime>
#include <algorithm>

#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"
#include <time.h>
#include <iomanip>
#include <unistd.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>

void* ProcessThread(void *_p);
void* RefineThread(void *_p);

void* thread(void *_p);

//bool lessThanSpec(const mspectrum &_l,const mspectrum &_r);
bool lessThanSpec(const mspectrum &_l,const mspectrum &_r)
{
	return _l.m_dMH < _r.m_dMH;
}
bool lessThanSequenceDebugUid(const msequence &_l,const msequence &_r)
{
	return _l.m_tUid < _r.m_tUid;
}

/*----------------------------------------------------------*/
mprocess **pProcess = NULL;
vector<msequenceCollection> mseqCol;
atomic<unsigned long> batchVecSize(0);
atomic<unsigned long> serilize_batchs(0);
bool reuse_serilize = true;

int pProCnt = 0;// count mprocess numbers
pthread_mutex_t lock_statistics = PTHREAD_MUTEX_INITIALIZER; 
pthread_mutex_t mtx_seqmap = PTHREAD_MUTEX_INITIALIZER; 
pthread_mutex_t lock_crc = PTHREAD_MUTEX_INITIALIZER; 
pthread_mutex_t lock_dup = PTHREAD_MUTEX_INITIALIZER; 

bool process_end = false;
size_t ProteinGlobalId = 0;
float crc_timeCnt = 0;

pthread_barrier_t bar_mpro;
pthread_barrier_t barrier_parallel;
pthread_barrier_t barrier_divide_wl;

int Protein_threads = 1;
int Spectra_divisions = 1;
int pipeline_rows = 0;
int pipeline_cols = 0;

pthread_t *tid = NULL;
class sequenceControl{
public:
	pthread_mutex_t mutex;
	pthread_cond_t cond;
	int   flag;
	float  cost;
	sequenceControl(void){
		flag = -1;
		pthread_mutex_init(&mutex,NULL);
		pthread_cond_init(&cond,NULL);
		cost = 0.0;
	}
	~sequenceControl(void){
		pthread_cond_destroy(&cond);
		pthread_mutex_destroy(&mutex);
	}
};
sequenceControl *statusBuffer = NULL;

class overheadDetail {
public:
	int seq_start;
	int seq_end;
	size_t workload;
	double thread_batch_overhead;
	double thread_parallel_overhead;
	double thread_pipeline_overhead;
	vector<double> block_overhead;
	vector<double> multiplex_overhead;
	overheadDetail(void){
		seq_start = 0;
		seq_end = 0;
		workload = 0;
		thread_batch_overhead = 0.0;
		thread_parallel_overhead = 0.0;
		thread_pipeline_overhead = 0;
	}
	~overheadDetail(void) {}
};

overheadDetail *overheadTable = NULL;
divideBatch *pBatchWL;

orderEle pSerializesArray[1<<16];
set<int> checkArraySet;

double thread_keep[128];
double thread_record[128];

/*----------------------------------------------------------*/
inline int set_cpu(int i){
	cpu_set_t mask;
	CPU_ZERO(&mask);

	CPU_SET(i,&mask);
	if(-1 == pthread_setaffinity_np(pthread_self(),sizeof(mask),&mask))
		return -1;
	return 0;
}


int main(int argc, char* argv[])
{
	/*
	* Check the argv array for at least one parameter.
	* mprocess checks the validity of the file.
	*/
	if(argc < 2 || argc > 1 && strstr(argv[1],"-L") == argv[1] || argc > 1 && strstr(argv[1],"-h") == argv[1])	{
		cout << "\n\nUSAGE: tandem filename\n\nwhere filename is any valid path to an XML input file.\n\n+-+-+-+-+-+-+\n";
		cout << "\nX! TANDEM " << VERSION << "\n";
		cout << "\nCopyright (C) 2003-2014 Ronald C Beavis, all rights reserved\n";
		cout << "This software is a component of the GPM  project.\n";
		cout << "Use of this software governed by the Artistic license.\n";
		cout << "If you do not have this license, you can get a copy at\n";
		cout << "http://www.perl.com/pub/a/language/misc/Artistic.html\n";
		cout << "\n+-+-+-+-+-+-+\n\npress <Enter> to continue ...";
		char *pValue = new char[128];
		cin.getline(pValue,127);
		delete pValue;
		return -1;
	}
	cout << "\nX! TANDEM " << VERSION << "\n\n";
	/*
	* Create an mprocess object array
	*/
	unsigned long lMaxThreads = 256;
#ifdef PROTEIN_GROUP
	pProcess = new mprocess*[lMaxThreads];
	mseqCol.reserve(1000000);
#else
	mprocess **pProcess = new mprocess*[lMaxThreads];//
#endif
	if(pProcess == NULL)	{
		cout << "An error was detected creating the processing objects.\nPlease contact a GPM administrator.\n";
		return -2;
	}

	pthread_t pThreads[lMaxThreads];

	unsigned long a = 0;
	while(a < lMaxThreads)	{
		pProcess[a] = NULL;
//		pHandle[a] = NULL;
//		pId[a] = NULL;
		a++;
	}
	pProcess[0] = new mprocess;
#ifdef PROTEIN_GROUP
	pProCnt++;
#endif
	cout << "Loading spectra" << endl;

	/*
	* Initialize the first mprocess object with the input file name.
	*/
	struct timeval load_start,load_end;
	gettimeofday(&load_start,NULL);
	char *pS = new char[1024];
	strcpy(pS,argv[1]);
	if(!pProcess[0]->load(pS))	{
		cout << "\n\nAn error was detected while loading the input parameters.\nPlease follow the advice above or contact a GPM administrator to help you.";
		delete pProcess[0];
		delete pProcess;
		return -4;
	}
	cout << " loaded.\n";
	if(pProcess[0]->m_vSpectra.size() == 0){
		cout << "No input spectra met the acceptance criteria.\n";
		cout.flush();
		delete pProcess[0];
		delete pProcess;
		return 1;
	}
#ifdef PLUGGABLE_SCORING
	cout << "Pluggable scoring enabled." << endl;
#endif
      /*
       * Start the mprocess object and wait for it to return.
       */
	unsigned long lThread = pProcess[0]->get_thread();
	unsigned long lThreads = pProcess[0]->get_threads();
	printf("thread:%d, %d",lThread,lThreads);
	if(lThreads	> lMaxThreads){
		lThreads = lMaxThreads;
	}
	if(pProcess[0]->m_vSpectra.size() <	lThreads){
		lThreads = (unsigned long)pProcess[0]->m_vSpectra.size();
		if(lThreads	< 1)	{
			lThreads = 1;
		}
		pProcess[0]->set_threads(lThreads);
	}

	int dCount = lThreads - 1;

	long lSpectra =	lThreads + (long)pProcess[0]->m_vSpectra.size()/lThreads;
	bool bSpectra =	true;
	cout << "Starting threads ." << endl;

	size_t tCount = pProcess[0]->m_vSpectra.size();
	sort(pProcess[0]->m_vSpectra.begin(),pProcess[0]->m_vSpectra.end(),lessThanSpec);//template sort
	float fMax = (float)pProcess[0]->m_vSpectra.back().m_dMH;
	float fZ = pProcess[0]->m_vSpectra.back().m_fZ;
	pProcess[0]->m_fMaxMass = fMax;
	pProcess[0]->m_fMaxZ = fZ;
	dCount = 0;
	pProcess[0]->serialize();
	cout << "Spectra matching criteria = " << (unsigned long)pProcess[0]->m_vSpectra.size() << "\n";
	gettimeofday(&load_end,NULL);
	double load_duration = (load_end.tv_sec - load_start.tv_sec) + (load_end.tv_usec-load_start.tv_usec)/1000000.0;
	cout<<"load duration "<< load_duration<<" seconds"<<endl;

	struct timeval proc_start, proc_end;
	gettimeofday(&proc_start,NULL);

	string strValue, strKey;
	strKey = "threads";
	pProcess[0]->m_xmlValues.get(strKey,strValue);
	Protein_threads = atoi(strValue.c_str());

	string of_pro = "system-profile.csv";
	pProcess[0]->of_profile.open(of_pro,ofstream::out | ofstream::app);
	pProcess[0]->of_profile << Protein_threads << ", " << load_duration << ", ";

//	strKey = "spectrum, divisions";
//	pProcess[0]->m_xmlValues.get(strKey,strValue);
//	Spectra_divisions = atoi(strValue.c_str());
	strKey = "pipeline, rows";
	pProcess[0]->m_xmlValues.get(strKey,strValue);
	pipeline_rows = atoi(strValue.c_str());
	strKey = "pipeline, cols";
	pProcess[0]->m_xmlValues.get(strKey,strValue);
	pipeline_cols = atoi(strValue.c_str());

	pthread_barrier_init(&bar_mpro,NULL,Protein_threads+1);
	pthread_barrier_init(&barrier_parallel,NULL,Protein_threads+1);
	pthread_barrier_init(&barrier_divide_wl,NULL,Protein_threads+1);

	tid = new pthread_t[Protein_threads];

	statusBuffer  = new sequenceControl[Protein_threads];
//	overheadTable = new overheadDetail[Protein_threads];
//	pBatchWL = new divideBatch[Protein_threads];

	cout << "Started. Computing models:" << endl;

	for(int i = 0; i < Protein_threads; ++i){
		pProcess[pProCnt] = new mprocess;
		pProcess[pProCnt]->set_class(pS);
		if(pthread_create(&tid[i],NULL,thread,(void*)pProCnt) != 0)
			cerr << "thread create fail" << endl;
		pProCnt++;
	}
	gettimeofday(&proc_end,NULL);
	double thread_config = (proc_end.tv_sec-proc_start.tv_sec) + (proc_end.tv_usec-proc_start.tv_usec)/1000000.0;
	cout<<"thread config : " << thread_config << " seconds"<<endl;
	pProcess[0]->of_profile << thread_config << ", ";

	gettimeofday(&proc_start,NULL);
	ProcessThread((void*)pProcess[0]); // main thread
	delete pS;
	for(int i = 0; i < Protein_threads; ++i){
		pthread_join(tid[i],NULL);
	}
#if 0
	for(int i = 0; i < Protein_threads;++i){
		int lib_start = i * batchVecSize / Protein_threads;
		int lib_end = (i+1) * batchVecSize / Protein_threads;
		for(int j = lib_start;j < lib_end;++j){
			(*(pProcess[i+1]->m_svrSequences.m_pCol)) = mseqCol[j];
		//	pProcess[i+1]->p_hydro = &pProcess[i+1]->thread_hydro[j-lib_start];
			pProcess[i+1]->p_hydro = pSerializesArray[j].p_batch;
		//	(*(tpro->m_svrSequences.m_pCol)) = mseqCol[readId];
			pProcess[i+1]->reuse_each_sequence();
		}
	}
#endif

	gettimeofday(&proc_end,NULL);
	double process_duration = (proc_end.tv_sec - proc_start.tv_sec)+(proc_end.tv_usec - proc_start.tv_usec) / 1000000.0; 
	cout << "serialize and parallel spend , "<< process_duration << ", ptandem thread , "<< Protein_threads <<endl;
	pProcess[0]->of_profile << process_duration << ", " << endl;
	for(int i = 0; i < Protein_threads;++i){
//		pProcess[0]->of_profile << i+1 << "," << thread_keep[i] << "," << thread_record[i] << endl;
	}

	pProcess[0]->clean_sequences();
//	cout << "\tsequences modelled = "<< (long)(pProcess[0]->get_protein_count()/1000.0 + 0.5) << " ks\n";
	cout << "\tsequences modelled = "<< (long)(pProcess[0]->get_protein_count()) << endl;

	pProcess[0]->merge_spectra();
	cout << "************************************" << endl;
	cout << " m_vseqBest   size       : " << pProcess[0]->m_vseqBest.size() << endl;
	cout << "************************************" << endl;

#if 10
	a = 1;
	/*
	* merge the results into the first object
	*/
	while(a < dCount)	{
		pProcess[0]->merge_map(pProcess[a]->m_mapSequences);
		pProcess[0]->merge_spectra(pProcess[a]->m_vSpectra);
		a++;
	}
	a = 1;
	pProcess[0]->load_sequences();
	while(a < dCount)	{
		pProcess[a]->merge_map(pProcess[0]->m_mapSequences);
		pProcess[a]->m_vseqBest = pProcess[0]->m_vseqBest;
		pProcess[a]->m_mapCrc.clear();
		a++;
	}
	/*
	* Report the contents of the mprocess objects into an XML file as described
	* in the input file.
	*/
	cout << "Model refinement:\n";
	time_t refine_start = clock();
	cout.flush();

	dCount = 0;

	RefineThread((void*)pProcess[dCount]);
	dCount++;
	/*
	* Initialize more mprocess objects, if lThread is not 0xFFFFFFFF, which signifies default single
	* threaded operation.
	*/
	if(lThread != 0xFFFFFFFF)	{
		while(dCount < lThreads)	{
			pthread_create(&pThreads[dCount],NULL,RefineThread,(void*)pProcess[dCount]);
			dCount++;
		}
	}
	/*
	* wait until all of the mprocess objects return.
	*/
	a = 1;
	/*
	* merge the results into the first object
	*/
	if(dCount > 1)	{
		cout << "Merging results:\n";
		cout.flush();
	}
	while(a < dCount)	{
		if(a == 1)	{
			cout << "\tfrom " << a+1;
		}
		else	{
			cout << a+1;
		}
		cout.flush();
		if(!pProcess[0]->add_spectra(pProcess[a]->m_vSpectra))	{
			cout << "adding spectra failed.\n";
		}
		pProcess[0]->merge_statistics(pProcess[a]);
		pProcess[a]->clear();
		pProcess[a]->m_mapSequences.clear();
		delete pProcess[a];
		pProcess[a] = NULL;
		a++;
	}
	if(dCount > 1)	{
		cout << "\n\n";
		cout.flush();
	}
	cout.flush();
	time_t refine_finish = clock();
	double refine_duration = (double)(refine_finish - refine_start)/CLOCKS_PER_SEC;
	cout<<"refine duration "<<refine_duration<<" seconds"<<endl;
	cout << "Creating report:\n";
	time_t report_start = clock();
	cout.flush();
	pProcess[0]->report();
	size_t tValid = pProcess[0]->get_valid();
	size_t tUnique = pProcess[0]->get_unique();
	double dE = pProcess[0]->get_error_estimate();
	unsigned long lE = (unsigned long)(0.5+dE);
	unsigned long lEe = (unsigned long)(0.5 + sqrt(dE));
	if(lEe == 0)	{
		lEe = 1;
	}
	if(dE <= 0.0)	{
		dE = 1.0;
	}
	cout << "\nValid models = " << (unsigned long)tValid << "\n";
	if(tUnique > 0)	{
		cout << "Unique models = " << (unsigned long)tUnique << "\n";
		cout << "Estimated false positives = " << lE << " &#177; ";
		cout << lEe << "\n";
	}
	lE = pProcess[0]->get_reversed();
	if(lE != -1)	{
		cout << "False positive rate (reversed sequences) = " << lE << "\n";
	}

	time_t report_finish = clock();
	double report_duration = (double)(report_finish - report_start)/CLOCKS_PER_SEC;
	cout << "report duration "<<report_duration<<" seconds"<<endl;
	cout << "\n\n";

#if 10
	ofstream pout("Protein_cnt.csv");
	for(const auto _p : pProcess[0]->protein_cnt){
		pout << _p.first << ", " << _p.second << endl;
	}
#endif
	/*
	* Delete the mprocess objects and exit
	*/
	a = 0;
	while(a < 16)	{
		if(pProcess[a] != NULL)	{
//			delete pProcess[a];
		}
		a++;
	}
#endif
	delete pProcess;
//	delete pId;
//	delete pHandle;
	
	return 0;

}
/*
 * Process thread is used to create the worker threads for each mprocess object
 */

void* ProcessThread(void *_p){
	((mprocess *)_p)->process();
	return (void*)0;
}

void* RefineThread(void *_p){
	((mprocess *)_p)->refine();
	return (void*)0;
}

void* thread(void *_p){

//	mprocess *tmp = (mprocess *)_p;// just pointer copy
	long pro_ID = ((long )_p);
	mprocess *tmp = pProcess[0];
	mprocess *tpro = pProcess[pro_ID];
	msequenceCollection item;
	bool endFlag = true;
	//synchronization
	pthread_barrier_wait(&bar_mpro);
	tpro->serialize();
	tpro->m_vSpectra = tmp->m_vSpectra;
	tpro->initialize();
	tpro->m_pScore->set_mtA_filter(0,tpro->m_vSpectra.size());
	int thread_i = pro_ID - 1;
	set_cpu(pro_ID);

	int lib_start = thread_i * batchVecSize / Protein_threads;
	int lib_end   = (thread_i+1) * batchVecSize / Protein_threads;

	struct timeval keep_start, keep_end;
	double keep_time;
	gettimeofday(&keep_start,NULL);
	int readId = 0;
	size_t peptideCount = 0;
	for(readId = lib_start;readId < lib_end;++readId){
		(*(tpro->m_svrSequences.m_pCol)) = mseqCol[readId];
		tpro->p_hydro = new vector<vector<hydrolyze_seq>>;
		tpro->hydro_each_sequence();

		tpro->thread_hydro.push_back(*(tpro->p_hydro));
		pSerializesArray[readId].p_batch = tpro->p_hydro;
		pSerializesArray[readId].workload = tpro->m_tPeptideCount;
		peptideCount += tpro->m_tPeptideCount;
		tpro->m_tPeptideCount = 0;
	}
#if 0
	while(1){
		endFlag = tmp->batchTaskQueue.Remove(readId);
		if(!endFlag)
			break;
		(*(tpro->m_svrSequences.m_pCol)) = mseqCol[readId];
		tpro->p_hydro = new vector<vector<hydrolyze_seq>>;
		tpro->hydro_each_sequence();

	/*	pthread_mutex_lock(&lock_crc);
		if(checkArraySet.find(readId) != checkArraySet.end())
			cout<<"array index : "<<readId<< endl;
		else
			checkArraySet.insert(readId);
		pthread_mutex_unlock(&lock_crc);*/

		pSerializesArray[readId].p_batch = tpro->p_hydro;
		pSerializesArray[readId].workload = tpro->m_tPeptideCount;
		peptideCount += tpro->m_tPeptideCount;
		tpro->m_tPeptideCount = 0;
	}
#endif
	__sync_fetch_and_add(&pProcess[0]->m_tPeptideScoredCount,tpro->m_tPeptideScoredCount);
	gettimeofday(&keep_end,NULL);
	keep_time = (keep_end.tv_sec-keep_start.tv_sec) + (keep_end.tv_usec-keep_start.tv_usec)/1000000.0;
	thread_keep[thread_i]=keep_time;

	pthread_barrier_wait(&barrier_parallel);
	pthread_barrier_wait(&barrier_divide_wl);

	gettimeofday(&keep_start,NULL);
#if 0
	while(1)
	{
		endFlag = tmp->batchTaskQueue.Remove(readId);
		if(!endFlag)
			break;
		tpro->p_hydro = pSerializesArray[readId].p_batch;
		(*(tpro->m_svrSequences.m_pCol)) = mseqCol[readId];
		tpro->parallel_reuse_each_sequence();
	}
#endif
#if 10
	if(thread_i >= pipeline_rows)
		return NULL;
	int Spectra_size = tpro->m_vSpectra.size();
	int libp_start = thread_i * batchVecSize / pipeline_rows;
	int libp_end   = (thread_i+1) * batchVecSize / pipeline_rows;

	for(int spec_i = 0;spec_i < pipeline_cols;++spec_i){
		int spec_start = (spec_i*Spectra_size) / pipeline_cols;
		int spec_end   = ((spec_i+1)*Spectra_size) / pipeline_cols;
		if(thread_i != 0){
			pthread_mutex_lock(&statusBuffer[thread_i-1].mutex);
			while(statusBuffer[thread_i-1].flag < spec_i){
				pthread_cond_wait(&(statusBuffer[thread_i-1].cond),&(statusBuffer[thread_i-1].mutex));
			}
			pthread_mutex_unlock(&statusBuffer[thread_i-1].mutex);
		}

		tpro->m_pScore->set_mtA_filter(spec_start,spec_end);
		for(int libp_i = libp_start;libp_i < libp_end;++libp_i){
			tpro->p_hydro = pSerializesArray[libp_i].p_batch;
			(*(tpro->m_svrSequences.m_pCol)) = mseqCol[libp_i];
			tpro->reuse_each_sequence();
		}

		if(thread_i < (pipeline_rows-1)){
			pthread_mutex_lock(&statusBuffer[thread_i].mutex);
			statusBuffer[thread_i].flag = spec_i;
			pthread_cond_signal(&statusBuffer[thread_i].cond);
			pthread_mutex_unlock(&statusBuffer[thread_i].mutex);
		}
	}
#endif

	gettimeofday(&keep_end,NULL);
	keep_time = (keep_end.tv_sec-keep_start.tv_sec) + (keep_end.tv_usec-keep_start.tv_usec)/1000000.0;
	thread_record[thread_i]=keep_time;
	return (void*)0;
}
//	size_t batch_start = serilize_batchs + ((pro_ID - 1) * batchVecSize) / Protein_threads;
//	size_t batch_end   = serilize_batchs + (pro_ID * batchVecSize) / Protein_threads ;
//	for(size_t batch_i = batch_start;batch_i < batch_end;++batch_i){
//		tpro->p_hydro = pSerializesArray[batch_i].p_batch;
//		(*(tpro->m_svrSequences.m_pCol)) = mseqCol[batch_i];
//		tpro->parallel_reuse_each_sequence();
//	}
