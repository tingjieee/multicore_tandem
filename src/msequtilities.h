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

#ifndef MSEQUENCEUTILITIES_H
#define MSEQUENCEUTILITIES_H

// File version: 2004-06-01

/*
 * msequtilities contains a set of constants and arrays for scoring and setting the molecular
 * mass of a peptide sequence. Biemann notation is used for fragment ions.
 */

#include "masscalc.h"
#include "mmotif.h"

typedef map<size_t,float> nMap;

class msequtilities
{
public:
	msequtilities(masscalc::massType _t);
	virtual ~msequtilities(void);

	masscalc m_calc;

	bool m_bPotential; // true if potential modifications are to be used
	bool m_bComplete; // true if complete modifications are to be used
	bool m_bPhosphoSerine;
	bool m_bPhosphoThreonine;
	bool m_bPhosphoTyrosine;
	bool m_bForcedMods;
	double m_dAmmonia; // mass of ammonia
	double m_dProton; // mass of proton
	double m_dWater; // mass of water
	double m_dHydrogen;

	float m_fNT; // mass to modify the N-terminal amino acid in a protein
	float m_fCT; // mass to modify the C-terminal amino acid in a protein
	double m_dA; // mass to add to a peptide mass to obtain an a ion mass
	double m_dB; // mass to add to a peptide mass to obtain a b ion mass
	double m_dC; // mass to add to a peptide mass to obtain a c ion mass
	double m_dX; // mass to add to a peptide mass to obtain an x ion mass
	double m_dY; // mass to add to a peptide mass to obtain a y ion mass
	double m_dZ; // mass to add to a peptide mass to obtain a z ion mass
	double m_dCleaveN; // mass added to the N-terminal residue by a cleavage reaction
	double m_dCleaveC; // mass added to the C-terminal residue by a cleavage reaction
	double m_dCleaveNdefault; // 1 (hydrolysis)
	double m_dCleaveCdefault; // 17 (hydrolysis)

	float *m_pfAaMass; // array of residue masses (less water), indexed by residue character code
	double *m_pdAaMass; // array of residue masses (less water), indexed by residue character code
	double *m_pdAaMod; // array of residue modification masses, indexed by residue character code
	double *m_pdAaPrompt; // array of residue modification masses, indexed by residue character code
	double *m_pdAaFullMod; // array of residue modification masses, indexed by residue character code
	float *m_pfAScore; // array of residue scores for a ions, indexed by residue character code
	float *m_pfBScore; // array of residue scores for b ions, indexed by residue character code
	float *m_pfYScore; // array of residue scores for y ions, indexed by residue character code
	float *m_pfXScore; // array of residue scores for x ions, indexed by residue character code
	float *m_pfA17Score; // array of residue scores for a-17 ions, indexed by residue character code
	float *m_pfB17Score; // array of residue scores for b-17 ions, indexed by residue character code
	float *m_pfY17Score; // array of residue scores for y-17 ions, indexed by residue character code
	float *m_pfX17Score;  // array of residue scores for x-17 ions, indexed by residue character code
	float *m_pfA18Score; // array of residue scores for a-18 ions, indexed by residue character code
	float *m_pfB18Score; // array of residue scores for b-18 ions, indexed by residue character code
	float *m_pfY18Score; // array of residue scores for y-18 ions, indexed by residue character code
	float *m_pfX18Score; // array of residue scores for x-18 ions, indexed by residue character code
	int *m_piAaDepth; // array of modification depth, indexed by residue character code
	int *m_piCtMod; // array of modification depth, indexed by residue character code
	int m_iSTYlimit;
	bool add_mod(const char _c,const size_t _v);
	bool clear_motifs(const bool _b);
	bool set_motifs(void);
	bool modify_all(string &_s);
	bool modify_c(const float _f);
	bool modify_maybe(string &_s);
	bool modify_annotation(string &_s);
	bool modify_motif(const string &_s);
	bool motif_set(const msequence &_s);
	string m_strDefaultMaybe;

#ifdef PROTEIN_GROUP
	msequtilities(const msequtilities& _util){
		m_bPotential = _util.m_bPotential ; // true if potential modifications are to be used
		m_bComplete = _util.m_bComplete; // true if complete modifications are to be used
		m_bPhosphoSerine = _util.m_bPhosphoSerine;
		m_bPhosphoThreonine = _util.m_bPhosphoThreonine;
		m_bPhosphoTyrosine = _util.m_bPhosphoTyrosine;
		m_bForcedMods = _util.m_bForcedMods;
		m_dAmmonia = _util.m_dAmmonia; // mass of ammonia
		m_dProton = _util.m_dProton; // mass of proton
		m_dWater = _util.m_dWater; // mass of water
		m_dHydrogen = _util.m_dHydrogen;

		m_fNT = _util.m_fNT; // mass to modify the N-terminal amino acid in a protein
		m_fCT = _util.m_fCT; // mass to modify the C-terminal amino acid in a protein
		m_dA = _util.m_dA ; // mass to add to a peptide mass to obtain an a ion mass
		m_dB = _util.m_dB; // mass to add to a peptide mass to obtain a b ion mass
		m_dC = _util.m_dC; // mass to add to a peptide mass to obtain a c ion mass
		m_dX = _util.m_dX; // mass to add to a peptide mass to obtain an x ion mass
		m_dY = _util.m_dY; // mass to add to a peptide mass to obtain a y ion mass
		m_dZ = _util.m_dZ; // mass to add to a peptide mass to obtain a z ion mass
		m_dCleaveN = _util.m_dCleaveN; // mass added to the N-terminal residue by a cleavage reaction
		m_dCleaveC = _util.m_dCleaveC; // mass added to the C-terminal residue by a cleavage reaction
		m_dCleaveNdefault = _util.m_dCleaveNdefault; // 1 (hydrolysis)
		m_dCleaveCdefault = _util.m_dCleaveCdefault; // 17 (hydrolysis)
	
		memcpy(m_pfAaMass,_util.m_pfAaMass,128*sizeof(float)); // array of residue masses (less water)
		memcpy(m_pdAaMass,_util.m_pdAaMass,128*sizeof(double)); // array of residue masses (less water)
		memcpy(m_pdAaMod,_util.m_pdAaMod,128*sizeof(double)); // array of residue modification masses
		memcpy(m_pdAaPrompt,_util.m_pdAaPrompt,128*sizeof(double)); // array of residue modification masses
		memcpy(m_pdAaFullMod,_util.m_pdAaFullMod,128*sizeof(double)); // array of residue modification masses
		
		memcpy(m_pfAScore,_util.m_pfAScore,128*sizeof(float)); // array of residue scores for a ions
		memcpy(m_pfBScore,_util.m_pfBScore,128*sizeof(float)); // array of residue scores for b ions
		memcpy(m_pfYScore,_util.m_pfYScore,128*sizeof(float)); // array of residue scores for y ions
		memcpy(m_pfXScore,_util.m_pfXScore,128*sizeof(float)); // array of residue scores for x ions
		
		memcpy(m_pfA17Score,_util.m_pfA17Score,128*sizeof(float)); // array of residue scores for a-17 ions
		memcpy(m_pfB17Score,_util.m_pfB17Score,128*sizeof(float)); // array of residue scores for b-17 ions
		memcpy(m_pfY17Score,_util.m_pfY17Score,128*sizeof(float)); // array of residue scores for y-17 ions
		memcpy(m_pfX17Score,_util.m_pfX17Score,128*sizeof(float)); // array of residue scores for x-17 ions
		
		memcpy(m_pfA18Score,_util.m_pfA18Score,128*sizeof(float)); // array of residue scores for a-18 ions
		memcpy(m_pfB18Score,_util.m_pfB18Score,128*sizeof(float)); // array of residue scores for b-18 ions
		memcpy(m_pfY18Score,_util.m_pfY18Score,128*sizeof(float)); // array of residue scores for y-18 ions
		memcpy(m_pfX18Score,_util.m_pfX18Score,128*sizeof(float)); // array of residue scores for x-18 ions

		memcpy(m_piAaDepth,_util.m_piAaDepth,128*sizeof(int)); // array of modification depth
		memcpy(m_piCtMod,_util.m_piCtMod,128*sizeof(int)); // array of modification depth
	
		m_iSTYlimit = _util.m_iSTYlimit;
		
		m_strDefaultMaybe = _util.m_strDefaultMaybe;
	
		m_vMotifs.clear();
		size_t a = 0;
		while( a < _util.m_vMotifs.size()){
			m_vMotifs.push_back(m_vMotifs[a]);
			a++;
		}
		m_mapMotifs.clear();
		if(_util.m_mapMotifs.size() > 0){
			m_mapMotifs = _util.m_mapMotifs;
		}
		m_mapMotifMods.clear();
		if(_util.m_mapMotifMods.size() > 0){
			m_mapMotifMods = _util.m_mapMotifMods;
		}
		m_mapMods.clear();
		if(_util.m_mapMods.size() > 0){
			m_mapMods = _util.m_mapMods;
		}
		m_mapNeutralLoss.clear();
		if(_util. m_mapNeutralLoss.size() > 0){
			m_mapNeutralLoss = _util.m_mapNeutralLoss;
		}
		m_bPotentialMotif = _util.m_bPotentialMotif;
		m_bSequenceMods = _util.m_bSequenceMods;
		m_bPrompt = _util.m_bPrompt;

		m_bIsModified = _util.m_bIsModified;
	}


	msequtilities& operator=(const msequtilities& _util){
		m_bPotential = _util.m_bPotential ; // true if potential modifications are to be used
		m_bComplete = _util.m_bComplete; // true if complete modifications are to be used
		m_bPhosphoSerine = _util.m_bPhosphoSerine;
		m_bPhosphoThreonine = _util.m_bPhosphoThreonine;
		m_bPhosphoTyrosine = _util.m_bPhosphoTyrosine;
		m_bForcedMods = _util.m_bForcedMods;
		m_dAmmonia = _util.m_dAmmonia; // mass of ammonia
		m_dProton = _util.m_dProton; // mass of proton
		m_dWater = _util.m_dWater; // mass of water
		m_dHydrogen = _util.m_dHydrogen;

		m_fNT = _util.m_fNT; // mass to modify the N-terminal amino acid in a protein
		m_fCT = _util.m_fCT; // mass to modify the C-terminal amino acid in a protein
		m_dA = _util.m_dA ; // mass to add to a peptide mass to obtain an a ion mass
		m_dB = _util.m_dB; // mass to add to a peptide mass to obtain a b ion mass
		m_dC = _util.m_dC; // mass to add to a peptide mass to obtain a c ion mass
		m_dX = _util.m_dX; // mass to add to a peptide mass to obtain an x ion mass
		m_dY = _util.m_dY; // mass to add to a peptide mass to obtain a y ion mass
		m_dZ = _util.m_dZ; // mass to add to a peptide mass to obtain a z ion mass
		m_dCleaveN = _util.m_dCleaveN; // mass added to the N-terminal residue by a cleavage reaction
		m_dCleaveC = _util.m_dCleaveC; // mass added to the C-terminal residue by a cleavage reaction
		m_dCleaveNdefault = _util.m_dCleaveNdefault; // 1 (hydrolysis)
		m_dCleaveCdefault = _util.m_dCleaveCdefault; // 17 (hydrolysis)
	
		memcpy(m_pfAaMass,_util.m_pfAaMass,128*sizeof(float)); // array of residue masses (less water)
		memcpy(m_pdAaMass,_util.m_pdAaMass,128*sizeof(double)); // array of residue masses (less water)
		memcpy(m_pdAaMod,_util.m_pdAaMod,128*sizeof(double)); // array of residue modification masses
		memcpy(m_pdAaPrompt,_util.m_pdAaPrompt,128*sizeof(double)); // array of residue modification masses
		memcpy(m_pdAaFullMod,_util.m_pdAaFullMod,128*sizeof(double)); // array of residue modification masses
		
		memcpy(m_pfAScore,_util.m_pfAScore,128*sizeof(float)); // array of residue scores for a ions
		memcpy(m_pfBScore,_util.m_pfBScore,128*sizeof(float)); // array of residue scores for b ions
		memcpy(m_pfYScore,_util.m_pfYScore,128*sizeof(float)); // array of residue scores for y ions
		memcpy(m_pfXScore,_util.m_pfXScore,128*sizeof(float)); // array of residue scores for x ions
		
		memcpy(m_pfA17Score,_util.m_pfA17Score,128*sizeof(float)); // array of residue scores for a-17 ions
		memcpy(m_pfB17Score,_util.m_pfB17Score,128*sizeof(float)); // array of residue scores for b-17 ions
		memcpy(m_pfY17Score,_util.m_pfY17Score,128*sizeof(float)); // array of residue scores for y-17 ions
		memcpy(m_pfX17Score,_util.m_pfX17Score,128*sizeof(float)); // array of residue scores for x-17 ions
		
		memcpy(m_pfA18Score,_util.m_pfA18Score,128*sizeof(float)); // array of residue scores for a-18 ions
		memcpy(m_pfB18Score,_util.m_pfB18Score,128*sizeof(float)); // array of residue scores for b-18 ions
		memcpy(m_pfY18Score,_util.m_pfY18Score,128*sizeof(float)); // array of residue scores for y-18 ions
		memcpy(m_pfX18Score,_util.m_pfX18Score,128*sizeof(float)); // array of residue scores for x-18 ions

		memcpy(m_piAaDepth,_util.m_piAaDepth,128*sizeof(int)); // array of modification depth
		memcpy(m_piCtMod,_util.m_piCtMod,128*sizeof(int)); // array of modification depth
	
		m_iSTYlimit = _util.m_iSTYlimit;
		
		m_strDefaultMaybe = _util.m_strDefaultMaybe;
		
		m_vMotifs.clear();
		size_t a = 0;
		while( a < _util.m_vMotifs.size()){
			m_vMotifs.push_back(m_vMotifs[a]);
			a++;
		}
		m_mapMotifs.clear();
		if(_util.m_mapMotifs.size() > 0){
			m_mapMotifs = _util.m_mapMotifs;
		}
		m_mapMotifMods.clear();
		if(_util.m_mapMotifMods.size() > 0){
			m_mapMotifMods = _util.m_mapMotifMods;
		}
		m_mapMods.clear();
		if(_util.m_mapMods.size() > 0){
			m_mapMods = _util.m_mapMods;
		}
		m_mapNeutralLoss.clear();
		if(_util. m_mapNeutralLoss.size() > 0){
			m_mapNeutralLoss = _util.m_mapNeutralLoss;
		}
		m_bPotentialMotif = _util.m_bPotentialMotif;
		m_bSequenceMods = _util.m_bSequenceMods;
		m_bPrompt = _util.m_bPrompt;

		m_bIsModified = _util.m_bIsModified;

		return *this;
	}
#endif

#ifdef PEPTIDE_GROUP
	bool pull_peptide(const msequtilities _util){
	return true;
	}

	bool push_peptide(msequtilities &_util){
return true;
	}

#endif

	double setPotentialMotif(const bool _b)
	{
		m_bPotentialMotif = _b;
		return m_bPotentialMotif;
	}


	double mztomh(double mz, float z)
	{
		return (mz - m_dProton)*z + m_dProton;
	}

	double getAaMass(char c, unsigned long t)
	{
		size_t tC = c;
		double dValue = m_pdAaMass[tC];
		dValue += m_pdAaMod[tC];
		dValue += m_pdAaFullMod[tC];
		dValue += m_pdAaPrompt[tC];

		if (m_bSequenceMods)
		{
			SMap::iterator itSeq = m_mapMods.find(t);
			if(itSeq != m_mapMods.end())
				dValue += itSeq->second;
		}

		return dValue;
	}
	int myPow(int x, int p)
	{
		if (p == 0) return 1;
		if (p == 1) return x;

		int tmp = myPow(x, p / 2);
		if (p % 2 == 0) return tmp * tmp;
		else return x * tmp * tmp;
	}
	vector<mmotif> m_vMotifs;
	map<size_t,size_t> m_mapMotifs;
	map<char,size_t> m_mapMotifMods;
	SMap m_mapMods;
	nMap m_mapNeutralLoss;
	bool m_bPotentialMotif;
	bool m_bSequenceMods;
	bool m_bPrompt;
	bool modify_n(const float _f);
	bool synthesis(const bool _b);
	bool set_aa_file(string &_p);
	bool is_modified();
	bool set_modified(const bool _b);

private:
	bool m_bIsModified;
	bool set_aa();
	bool set_a(const char _c,const float _f);
	bool set_b(const char _c,const float _f);
	bool set_x(const char _c,const float _f);
	bool set_y(const char _c,const float _f);
};
#endif
