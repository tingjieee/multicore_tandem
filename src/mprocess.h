/*
 Copyright (C) 2003 Ronald C Beavis, all rights reserved
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

#ifndef MPROCESS_H
#define MPROCESS_H

// File version: 2003-08-01
// File version: 2004-03-01
// File version: 2004-11-01

typedef map<size_t,string> SEQMAP;

#include <sys/timeb.h>
#include <ctime>
#include "mcleave.h"
#include "mspectrumcondition.h"
#include "mreport.h"
#include "mzid_report.h"
#include "msequenceserver.h"
#include "msequencecollection.h"
#include "mplugin.h"
#include "msemistate.h"
#include "mscore.h"
#include "mrefine.h" 
#include <set>
#include <unordered_map>
#include <semaphore.h>
#include <pthread.h>
#include <sys/time.h>
#include <list>
#include <atomic>
#include <mutex>
#include <bitset>
typedef map<uint64_t,vector<msequence> > MSEQMAP;

/*
 * the process object coordinates the function of tandem. it contains the information
 * loaded from the input XML file in the m_xmlValues object and performance
 * information in the m_xmlPerformance object. The mass spectra to be analyzed are
 * in the m_vSpectra vector container. A set of input parameters are used to
 * initialize constants that are used in processing the mass spectra.
 */
class merrors
{
public:
	merrors(void) { 
		m_bPpm = true; 
		m_bIsotope = false;
		m_dPlus = 20.0; 
		m_dMinus = 20.0;
	}
	virtual ~merrors(void) { }
	bool m_bPpm;
	bool m_bIsotope;
	double m_dPlus;
	double m_dMinus;
	bool check(double _s,double _m)	{
		double dDelta = _s - _m;
		double dPlus = m_dPlus;
		double dMinus = m_dMinus;
		if(m_bPpm)	{
			dPlus *= _m*1.0e-6;
			dMinus*= _m*1.0e-6;
		}
		if(dDelta < 0.0)	{
			if(dDelta >= dMinus)	{
				return true;
			}
		}
		else	{
			if(dDelta <= dPlus)	{
				return true;
			}
		}
		if(!m_bIsotope)	{
			return false;
		}
		if(_s > 1000.0)	{
			dDelta -= 1.008664916;
			if(dDelta < 0.0)	{
				if(dDelta >= dMinus)	{
					return true;
				}
			}
			else	{
				if(dDelta <= dPlus)	{
					return true;
				}
			}
		}
		if(_s > 1500.0)	{
			dDelta -= 1.008664916;
			if(dDelta < 0.0)	{
				if(dDelta >= dMinus)	{
					return true;
				}
			}
			else	{
				if(dDelta <= dPlus)	{
					return true;
				}
			}
		}
		return false;
	}
};

class mprocesslog
{
public:
	mprocesslog() { }
	virtual ~mprocesslog() { }
	bool open(string &_s)	{
		m_ofLog.open(_s.c_str());
		if(m_ofLog.fail())	{
			return false;
		}
		return true;
	}
	bool log(string &_m)	{
		if(!m_ofLog.is_open())	{
			return false;
		}
		time_t tValue;
		time(&tValue);
		struct tm *tmValue = localtime(&tValue);
		char pLine[256];	
		strftime(pLine, 255,"%Y-%m-%d %H:%M:%S",tmValue);
		m_ofLog << pLine << "\t" << _m.c_str() << "\n";
		m_ofLog.flush();
		return true;
	}
	bool log(const char *_m)	{
		if(!m_ofLog.is_open())	{
			return false;
		}
		string strValue = _m;
		return log(strValue);
	}
	bool close()	{
		if(!m_ofLog.is_open())	{
			return false;
		}
		m_ofLog.close();
		return true;
	}
private:
	ofstream m_ofLog;
};

class item {
public:
	std::atomic<int> flag; // 0->writeable, 1->readable
	std::atomic<int> vectorindex;
	std::mutex mtx;
	int pad01;
	int pad02;
	int pad03;
	int pad04;
	int pad05;
	int pad06;
	int pad07;

	item(void){
		//flag = 0;
		flag.store(0);
		vectorindex.store(0);
	}
//	item(int i,int long ii) : flag(i),vectorindex(ii){}
	~item(void){
	}
};

class buf {
public:
	unsigned long pad01;
	unsigned long pad02;
	unsigned long pad03;
	unsigned long pad04;
	unsigned long pad05;
	unsigned long pad06;
	unsigned long pad07;
	unsigned long tail;
	unsigned long pad11;
	unsigned long pad12;
	unsigned long pad13;
	unsigned long pad14;
	unsigned long pad15;
	unsigned long pad16;
	unsigned long pad17;
	std::atomic<bool> produceStop;
	//bool produceStop;
	unsigned long pad1;
	unsigned long pad2;
	unsigned long pad3;
	unsigned long pad4;
	unsigned long pad5;
	unsigned long pad6;
	unsigned long pad7;
	size_t bufSize;
	item *p_item;
	unsigned long pad21;
	unsigned long pad22;
	unsigned long pad23;
	unsigned long pad24;
	unsigned long pad25;
	unsigned long pad26;
	unsigned long pad27;
	//float buf_Inserttime_cnt;
	//long buf_Removetime_cnt;

	buf(void)
	{
		bufSize = 128;
		p_item = new item[bufSize];
		tail = 0;
		produceStop.store(false);
	//	buf_Inserttime_cnt = 0;
	//	buf_Removetime_cnt = 0;
	}
	~buf(void){
		delete [] p_item;
	}

	bool Insert(const unsigned long _id)
	{
		unsigned long currentWriteIndex = 0;
		int flag= 0;
		do {
			currentWriteIndex = tail;
			//flag = p_item[currentWriteIndex % (bufSize)].flag.load();
			flag = p_item[(currentWriteIndex & (bufSize-1))].flag.load();
			if(flag == 0)
				break;
			++tail;
		} while(1);
		//p_item[currentWriteIndex & (bufSize-1)].vectorindex = _id;
		p_item[currentWriteIndex & (bufSize-1)].vectorindex.store(_id);
		//p_item[currentWriteIndex & (bufSize-1)].flag = 1;
		p_item[currentWriteIndex & (bufSize-1)].flag.store(1);

		return true;
	}
	bool Remove(int& readableIndex){
		int flag = 0;
		size_t i = 0;
		size_t circle = 0;
		int readFlag = 1, writeableFlag = 0;
		do{
			for(i = 0; i < bufSize;++i){
				p_item[i].mtx.lock();
				if(p_item[i].flag.load() == readFlag){
					readableIndex = p_item[i].vectorindex.load();
					p_item[i].flag.store(writeableFlag);
					p_item[i].mtx.unlock();
					return true;
				}
				p_item[i].mtx.unlock();
			}
			++circle;
			if(circle > 1 && produceStop.load())
				return false;
		}while(1);

		return true;
	}
#if 0
	//optimization next step , random the first index, avoid hot zone 
	bool Remove(int& readableIndex)
	{
		int flag = 0;
		size_t i = 0;
		size_t circle = 0;
		//timeval remove_start,remove_end;
		//long spend = 0;
		int readFlag = 1;
		do {
			i = 0;
		//	circle = 0;
			for( ; i < bufSize; ++i){
				if(flag == p_item[i].flag.load())
					break;
			//	if(flag == readFlag)	break;
			//	else	++circle;
			}
			++circle;
			//search all the ringbuffer,not found and done()
			//if(circle > bufSize && produceStop.load())
			//if(circle > 1 && produceStop.load())
			if(i >= bufSize && produceStop.load())
				return false;
		//} while(!__sync_bool_compare_and_swap(&(p_item[i].flag), 1, 0));
		}while(!(p_item[i].flag.compare_exchange_weak(readFlag,2)));
		//readableIndex = p_item[i].vectorindex;
		readableIndex = p_item[i].vectorindex.load();
		p_item[i].flag.store(0);
//	gettimeofday(&remove_end,NULL);
//	spend  =  (remove_end.tv_sec-remove_start.tv_sec)*1000000+(remove_end.tv_usec-remove_start.tv_usec);
//	__sync_fetch_and_add(&buf_Removetime_cnt,spend);

		return true;
	}
#endif
};



class unit {
public:
	long IonCount;
	float fScore;
	float fHyper;
	int chBCount;
	int chYCount;
//	std::array<unsigned long,16> plCount;
//	std::array<float,16> pfScore;
	unit(void) {
		IonCount=0;
		fScore = 0;
		fHyper = 0;
		chBCount=0;
		chYCount=0;
//		long a = 0;
//		while(a < 16){
//			plCount[a] = 0;
//			pfScore[a] = 0;
//			++a;
//		}
	}
	unit(const unit &_rhs){
		IonCount = _rhs.IonCount;
		fScore = _rhs.fScore;
		fHyper = _rhs.fHyper;
		chBCount=_rhs.chBCount;
		chYCount=_rhs.chYCount;
//		for(size_t i=0;i<16;++i){
//			plCount[i] = _rhs.plCount[i];
//			pfScore[i] = _rhs.pfScore[i];
//		}
	}
	unit& operator= (const unit &_rhs){
		if(&_rhs != this){
			IonCount = _rhs.IonCount;
			fScore = _rhs.fScore;
			fHyper = _rhs.fHyper;
			chBCount=_rhs.chBCount;
			chYCount=_rhs.chYCount;
//			for(size_t i=0;i<16;++i){
//				plCount[i] = _rhs.plCount[i];
//				pfScore[i] = _rhs.pfScore[i];
//			}
		}
		return *this;
	}
	~unit(void){
	}
};

class singleSpectrum {
public:
	double seqMH;
	long equals;
//	vector<theorSpectrum> tSpectrum;
	vector<size_t> m_plEqualsS;
	vector<unit> m_plRes;
//	vector<float> m_plHyper;
//	vector<float> m_plScore;
	singleSpectrum(void) {
		seqMH = 0;
		equals= 0;
	//	m_plEqualsS.reserve(1<<8);
	//	m_plRes.reserve(1<<8);
	//	m_plHyper.reserve(1<<8);
	//	m_plScore.reserve(1<<8);
	//	tSpectrum.clear();
	}
	~singleSpectrum(void){
	}
};

class onePeptide {
public:
	long start;
	long end;
	long missedCleaves;
//	long hydrStage;
	//multiset<double> seqMH;
	vector<singleSpectrum> Spec;
	onePeptide(void) {
		start = 0;
		end = 0;
		missedCleaves = 0;
//		hydrStage = 0;
	}
	~onePeptide(void) {
	}
};
class hydrolyze_seq {
public:
	long start_0;
	long end_0;
	long lastCleave;

	vector<onePeptide> seqStage_1;
	vector<onePeptide> seqStage_2;
	hydrolyze_seq(void)  { 
		start_0 = 0;
		end_0 = 0;
		lastCleave = 0;
	//	seqStage_1.reserve(1<<8);
	//	seqStage_2.reserve(1<<8);
	}
	~hydrolyze_seq(void) {
//		for(size_t i = 0; i < seqStage_1.size();)
	}
};

class orderEle {
public:
	vector<vector<hydrolyze_seq> > *p_batch;
	long workload;
	size_t pad01;
	size_t pad02;
	size_t pad03;
	size_t pad04;
	size_t pad05;
	size_t pad06;
	size_t pad07;
	size_t pad08;
	orderEle(void) {
		p_batch = NULL;
		workload = 0;
	}
	~orderEle(void){
	}
	orderEle& operator= (const orderEle& _mpl){
		if(&_mpl != this){
			p_batch = _mpl.p_batch;
			workload = _mpl.workload;
		}
		return *this;
	}
};
class divideBatch {
public:
	int start;
	int end;
	size_t workload;
	divideBatch(void){
		start = 0;
		end = 0;
		workload = 0;
	}
	~divideBatch(void){}
};

class mprocess
{
public:
	
#ifdef PROTEIN_GROUP
	long  crc_time_thread;
	double hydro_t;
	double multi_t;
	double recored_score_t;
	double score_t;
	double total_seq_time;
	double keep_time;
//	double create_score_t;
//	long workload;
//	vector<float> m_plHyper;
//	vector<float> m_plScore;
//	vector<size_t> m_plEqualsS;
	
	bool hydrolyze_repeat;
	map<size_t, size_t> protein_cnt;
	vector<vector<hydrolyze_seq>> m_hydro;
	vector<vector<vector<hydrolyze_seq>>> thread_hydro;
	vector<vector<hydrolyze_seq>> *p_hydro;
//	list<vector<hydrolyze_seq> > l_hydro;
//	vector<vector<hydrolyze_seq> > m_hydro_1;
	size_t s_rid;
	size_t s_wid;
//	size_t details_times;
	float Max_fZ;
	//vector<unsigned int> v_uiType;

	vector<orderEle> tmpOrder;

	buf batchTaskQueue;
	bitset<10000> spec_3k_vec;
	
	ofstream of_profile;

#endif

	mprocess(void);
	virtual ~mprocess(void);
	bool serialize(void); //serializes the m_vMI data elements in m_vSpectra onto the disk to save space during calculation
	bool restore(void); //restores the m_vMI data elements from disk
	bool removeMI(void);
	vector<string> m_vstrPaths;
	mprocesslog m_prcLog;
	bool add_spectra(vector<mspectrum> &_v); // adds the spectra contained in _v to m_vSpectra
	bool clear(); // clears a selection of vectors in the mprocess object
	unsigned long get_thread(); // retrieves the number of the current thread (0 - 15)
	unsigned long get_threads(); // retrieves the total number of threads running (1 - 16)
	size_t get_protein_count(); // gets the total number of proteins read
	size_t get_peptide_count(); // gets the total number of peptides used
	size_t get_total_residues(); // gets the total number of residues read
	size_t get_valid(); // gets the number of valid peptide models
	double get_error_estimate() {return m_dEsum;}
	size_t get_unique(); // gets the number of valid peptide models
	long get_reversed(); // gets the number of reversed peptide models
	double get_threshold(); // gets the number of valid peptide models
	bool load(const char *_f,mprocess *_p = NULL); // loads input parameters
	bool set_class(const char *_f,mprocess *_p = NULL); // loads input parameters
	bool load_saps(mprocess *_p); // loads sap information, if it exists
	bool load_annotation(mprocess *_p); // loads sequence annotation information, if it exists
	virtual bool merge_spectra(); // adds externally generated mspectrum vector to m_vSpectra
	virtual bool merge_spectra(vector<mspectrum> &_s); // adds externally generated mspectrum vector to m_vSpectra
	bool merge_map(SEQMAP &_s); // adds externally generated mspectrum vector to m_vSpectra
	bool merge_statistics(const mprocess *_p); // adds externally generated spectra to this object
	bool process(void); // performs identifications based on the input parameters
	bool initialize(void); // performs identifications based on the input parameters
	bool refine(void); // controls the protein model refinement process
	bool report(void); // produces the XML output report, using an mreport class
	bool pyro_check(const char _c);
	bool pyro_reset(void);
	bool score_each_sequence(); // generates a score for all sequences in the m_svrSequences object
	bool hydro_each_sequence(); // generates a hydrolyze for all sequences in the m_svrSequences object
	bool reuse_each_sequence(); // generates a hydrolyze for all sequences in the m_svrSequences object
	bool parallel_reuse_each_sequence(); // generates a hydrolyze for all sequences in the m_svrSequences object
	bool set_threads(const unsigned long _t); // sets the object's thread number
	bool set_thread(const unsigned long _t); // sets the object's thread number
	bool dynamic_parent_ion_selection(const double _d,const double _maxe); //selects a parent ion mass window for high resolution spectra dalton space
	bool dynamic_parent_ion_selection_ppm(const double _d, const double _maxe); //selects a parent ion mass window for high resolution spectra in ppm space
	int set_round(const int _r)	{ m_iCurrentRound = _r; return m_iCurrentRound;}
	virtual bool load_sequences();
	XmlParameter m_xmlPerformance; // stores process performance parameters
	XmlParameter m_xmlValues; // store process input parameters
	vector<mspectrum> m_vSpectra; // store spectra to be analyzed
	SEQMAP m_mapSequences; // a map containing all of the protein sequences discovered, indexed by their m_tUid value
	vector<msequence> m_vseqBest; // a vector of msequences used in the model refinement process
#ifdef DEBUG
	vector<msequence> m_debug;
#endif
	vector<string> m_vstrModifications; //a vector containing the strings defining fixed modifications for a protein
	size_t m_tRefineModels; // total number of models generated by refinement
	size_t m_tRefineInput; // total number of sequences included in a refinement session
	size_t m_tRefinePartial; // the number of models discovered to have partial cleavage
	size_t m_tRefineUnanticipated; // the number of models discovered to have unanticpated cleavage
	size_t m_tRefineNterminal; // the number of models discovered to have modified N-terminii
	size_t m_tRefineCterminal; // the number of models discovered to have modified C-terminii
	size_t m_tRefinePam; // the number of models discovered to have point mutations
	double m_dRefineTime; // the time required to perform a refinement
	size_t m_tActive;	// total number of models remaining after each refinement step
	bool m_bRefineCterm;  //true if processing 'refine, potential C-terminus modifications'. Set in mrefine::refine and 
						//checked in score(), so the start position can be set to the length of the protein sequence
						// minus the value for 'refine, potential N-terminus modification position limit' before performing cleavage
	vector<int> m_viQuality; // contains the data quality scoring vector
	//map<uint64_t,size_t> m_mapCrc;
	unordered_map<uint64_t,size_t> m_mapCrc;
	size_t m_tTemp;
	bool m_bReversedOnly;
	bool m_bSaps;
	bool m_bAnnotation;
	bool m_bMinimalAnnotation;
	bool m_bSerialize;
	bool m_bCheckNg;
	bool m_bSkyline;
	bool m_bUseCrc;
	string m_strSkyline;
	double m_dNt;
	double m_dNtAve;
	double m_dNg;
	double m_dNgAve;
	float m_fMaxMass;
	float m_fMaxZ;
	enum	{
		I_Y =	0x01,
		I_B =	0x02,
		I_X =	0x04,
		I_A =	0x08,
		I_C =	0x10,
		I_Z =	0x20,
	} ion_type; // enum for referencing information about specific ion types.

#ifdef PROTEIN_GROUP
	char* Get_mSeq(void) { return m_pSeq;}
	
//	virtual bool clean_sequences(void); // remove sequences no longer required from m_mapSequences
//	void *thread(void *vargp);

#endif


protected:
	/*
	* refinement classes are declared friend classes so they can access private member variables and functions. 
	* If a new refinement class is added to the project, include a declaration here
	*/
	friend class mrefine;
	friend class mpmods;
	friend class mxxcleavage;
	friend class mtermmods;
	friend class mpam;
#ifdef PROTEIN_GROUP
public:
#else
protected:
#endif
	string m_strLastMods;
	int m_iCurrentRound;
	bool m_bPermute;
	bool m_bPermuteHigh;
	bool m_bCrcCheck;
	bool m_bQuickAcetyl;
	bool m_bQuickPyro;
	double m_dEsum;
	bool m_bChargeDither;
	float m_fSavResolve;
	bool m_bSavResolve;
	set<size_t> m_setRound;
	vector<string> m_vstrSaps;
	vector<string> m_vstrMods;
	map<string,string> m_mapAnnotation;
	map<string,string>* m_pmapAnnotation;
	msemistate m_semiState; // maintains the state of the semi-enzymatic cleavage state machine
	mpyrostate m_pyroState; // maintains the state of the pyrolidone carboxylic acid detection state machine
	merrors m_errValues;
	double m_dSearchTime; // total time elapsed during a protein modeling session process
	long m_lIonCount; // minimum sum of detected ions that are significant enough to store a sequence 
	unsigned long m_lThread; // thread number of this object
	unsigned long m_lThreads; // the total number of threads current active
	long m_lReversed; // the total number of peptides found where the reversed sequence was better than the forward sequence
	double m_dThreshold; // the current expectation value threshold
	double m_tContrasted; // the number of spectra subtracted using contrast angle redundancy detection
	long m_lStartMax; // set the maximum distance from N-terminus for a peptide
					  // normally set at an impossibly large value = 100,000,000
					  // for ragged N-terminii with potential modifications, set at a low but plausible value = 50
	long m_lCStartMax;
	char *m_pSeq; // a character pointer, used for temporary sequence information
	bool m_bUn; // if true, cleave at all residues. if false, use cleavage specification in data input.
	bool m_bUseHomologManagement; // set to true to use homologue management 
	size_t m_tMinResidues; // the minimum peptide length that will be scored
	size_t m_tMaxResidues; // the maximum peptide length that will be scored
	size_t m_tMissedCleaves; // the maximum number of cleavage sites that can be missed
	size_t m_tPeptideCount; // the total number of peptide sequences generated during a process
	size_t m_tPeptideScoredCount; // the total number of peptide sequences scored during a process
	size_t m_tProteinCount; // the total number of protein sequences considered during a process
	size_t m_tSpectra; // the total number of spectra being modeled
	size_t m_tSpectraTotal; // the total number of spectra in the input file
	size_t m_tValid; // the number of valid peptide models
	size_t m_tTotalResidues; // the number of residues read
	size_t m_tSeqSize; // current length of the m_pSeq character array
	size_t m_tUnique; // the number of unique peptides found in a result
	string m_strOutputPath; // the path name of the XML output file
	mcleave m_Cleave; // the specification for a cleavage peptide bond
	msequence m_seqCurrent; // the msequence object that is currently being scored

#ifdef PROTEIN_GROUP
#ifdef X_P3
	p3msequenceServer m_svrSequences; // the msequenceServer object that provides msequences to msequenceCollection
#else
	msequenceServer m_svrSequences; // the msequenceServer object that provides msequences to msequenceCollection
#endif
#endif

	mspectrumcondition m_specCondition; // the mspectrumcondition object that cleans up and normalized
										// spectra for further processing
	mscore* m_pScore; // the object that is used to score sequences and spectra
	mrefine* m_pRefine; // the object that is used to refine models

	bool charge(void); // adds additional charge states to the list of spectra, resulting in +1, +2 and +3 all being represented
	virtual bool clean_sequences(void); // remove sequences no longer required from m_mapSequences
	double dot(const size_t _f,const size_t _s,const float _r,const bool _t); // calculated inner product between two spectrum vectors
	bool create_rollback(vector<mspectrum> &_v); // create a temporary rollback vector
	virtual bool create_score(const msequence &_s,const size_t _v,const size_t _w, const long _m,const bool _p); // generates scores for a sequence
	virtual bool create_keep_score(const msequence &_s,const size_t _v,const size_t _w, const long _m,const bool _p,vector<unit> &_plRes); // generates scores for a sequence
	bool multiplex_create_score(const msequence &_s,const size_t _v,const size_t _w, const long _m,const bool _p,const vector<unit> &_plRes); // generates scores for a sequence
	bool parallel_multiplex_create_score(const msequence &_s,const size_t _v,const size_t _w, const long _m,const bool _p,const vector<unit> &_plRes); // generates scores for a sequence
//	bool create_score(const msequence &_s,const size_t _v,const size_t _w, const long _m,const bool _p,const vector<theorSpectrum> &_Spec); // generates scores for a sequence
//	virtual bool create_hydro(vector<theorSpectrum> &seq_hydro); //generates hydro for a sequence
	double expect_protein(const unsigned long _c,const unsigned long _t,
							const unsigned long _n,const double _d); // assigns an expectation value for a protein
	bool load_best_vector(void); // creates a new m_vseqBest vector and marks assigned spectra inactive
	bool mark_repeats(void); // checks for repeated assignments of the same peptide sequence
	bool modify(void); // sets the initial residue modification parameters in the m_seqUtil object
	bool refine_model(void); // the method that refines protein models
	bool report_all(); // create a report that contains all peptide models, regardless of expectation value
	bool report_valid(const double _d); // create a report that only contains valid peptide models (expect < m_dThreshold)
	bool report_stochastic(const double _d); // create a report that contains only stochastic peptide models (expect > m_dThreshold)
	bool report_expect(const double _m); // calculates expectation values for the report
	bool report_sort(void); // sorts peptide models for use in the report
	bool residues(); // estimates the minimum number of residues that can possibly be contained in a peptide
	bool rollback(vector<mspectrum> &_v,const double _m,const double _f);
	bool score(const msequence &_s); // generates scores of an msequence object
	bool hydro(const msequence &_s); // generates hydro of an msequence object
	bool multiplex(const msequence &_s); // generates hydro of an msequence object
	bool parallel_multiplex(const msequence &_s); // generates hydro of an msequence object
	bool score_single(const msequence &_s); // generates scores of an msequence object
	bool hydro_single(const msequence &_s); // generates scores of an msequence object
	bool multiplex_single(const msequence &_s); // generates scores of an msequence object
	bool parallel_multiplex_single(const msequence &_s); // generates scores of an msequence object
	bool score_terminus(const string &_s); // attempts to find N or C terminal modified peptides
	bool score_terminus_single(const string &_s); // attempts to find N or C terminal modified peptides
	bool spectra(void); // loads the m_vSpectra object using the m_specCondition object
	bool spectra_force(string &_t,string &_v); // forces the spectrum loader to use a specified file type
	bool subtract(void); // remove redundant mass spectra
	uint64_t crc(const string &_s);
	uint64_t *m_pCrcTable;
	bool initialize_crc(void);
	MSEQMAP m_mapDups;
//	SEQMAP m_mapTest;	used to test the utility of the crc function, see commented code in mprocess::score
	size_t m_tDuplicates;
	size_t m_tDuplicateIds;
	bool insert_dups(void);
	virtual bool taxonomy(void); // loads the taxonomy setting into the m_svrSequences object
};

#include "p3mprocess.h"

#endif
