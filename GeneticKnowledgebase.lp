% Genetic Knowledgebase
% Contains some logical rules about nucleotide, amino acids, and proteins 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#script (python) 
import clingo
def stringLen(s):
    return len(s.string)
def letteri(thes,thei):
    s = thes.string
    i = thei.number  - 1 
    return "<RANGE ERR>" if i > len(s) else s[i]
def concatfun(s1,s2):
    return s1.string + s2.string
#end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%       Part 1       %%%

% define nucleotide bases
base(adenine).
base(thymine).
base(cytosine).
base(guanine).

% define base types
purine(adenine).
pyrimidine(thymine).
pyrimidine(cytosine).
purine(guanine).

% assert that if a base is a purine, then it is not a pyrimidine
-purine(BASE) :- pyrimidine(BASE).
-pyrimidine(BASE) :- purine(BASE).

% assert adenine bonds with thymine
bondsWith(adenine, thymine).

% define the general bondsWith/2 relation
bondsWith(BASE1, BASE2) :- 
	base(BASE1),
	base(BASE2),
	purine(BASE1), 
	pyrimidine(BASE2),
	pyrimidine(BASE3),
	not bondsWith(BASE1, BASE3),
	BASE2 != BASE3.

% declare symmetry of the bondsWith/2 relation
bondsWith(BASE1, BASE2) :-
	bondsWith(BASE2, BASE1).
	
% display answer
answer(WHAT):-bondsWith(cytosine,WHAT).
#show answer/1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%      Part Two      %%%

% codons are combinations of any 3 bases
% codon/1
codon(@concatfun(@concatfun(CODE1, CODE2), CODE3)) :-
    nucleotideCode(BASE1,CODE1),
    nucleotideCode(BASE2,CODE2),
    nucleotideCode(BASE3,CODE3),
    base(BASE1),
    base(BASE2),
    base(BASE3).

% codons encode amino acids or stops
% codon/2
codon(phenylalanine, "TTT").
codon(phenylalanine, "TTC").
codon(leucine, "TTA").
codon(leucine, "TTG").
codon(leucine, "CTT").
codon(leucine, "CTC").
codon(leucine, "CTA").
codon(leucine, "CTG").
codon(isoleucine, "ATT").
codon(isoleucine, "ATC").
codon(isoleucine, "ATA").
codon(methionine, "ATG").
codon(valine, "GTT").
codon(valine, "GTC").
codon(valine, "GTA").
codon(valine, "GTG").
codon(serine, "TCT").
codon(serine, "TCC").
codon(serine, "TCA").
codon(serine, "TCG").
codon(proline, "CCT").
codon(proline, "CCC").
codon(proline, "CCA").
codon(proline, "CCG").
codon(threonine, "ACT").
codon(threonine, "ACC").
codon(threonine, "ACA").
codon(threonine, "ACG").
codon(alanine, "GCT").
codon(alanine, "GCC").
codon(alanine, "GCA").
codon(alanine, "GCG").
codon(tyrosine, "TAT").
codon(tyrosine, "TAC").
codon(histidine, "CAT").
codon(histidine, "CAC").
codon(glutamine, "CAA").
codon(glutamine, "CAG").
codon(asparagine, "AAT").
codon(asparagine, "AAC").
codon(lysine, "AAA").
codon(lysine, "AAG").
codon(asparticAcid, "GAT").
codon(asparticAcid, "GAC").
codon(glutamicAcid, "GAA").
codon(glutamicAcid, "GAG").
codon(cysteine, "TGT").
codon(cysteine, "TGC").
codon(tryptophan, "TGG").
codon(arginine, "CGT").
codon(arginine, "CGC").
codon(arginine, "CGA").
codon(arginine, "CGG").
codon(arginine, "AGA").
codon(arginine, "AGG").
codon(serine, "AGT").
codon(serine, "AGC").
codon(glycine, "GGT").
codon(glycine, "GGC").
codon(glycine, "GGA").
codon(glycine, "GGG").

% identify start codons
startCodon("TTG").
startCodon("ATG").
startCodon("GTG").

% identify stop codons
stopCodon("TAA").
stopCodon("TAG").
stopCodon("TGA").

% nucleotide codes 
nucleotideCode(adenine,"A").
nucleotideCode(guanine,"G").
nucleotideCode(cytosine,"C").
nucleotideCode(thymine,"T").

% amino acid codes
aminoAcidCode(phenylalanine, "F").
aminoAcidCode(leucine, "L").
aminoAcidCode(isoleucine, "I").
aminoAcidCode(methionine, "M").
aminoAcidCode(valine, "V").
aminoAcidCode(serine, "S").
aminoAcidCode(proline, "P").
aminoAcidCode(threonine, "T").
aminoAcidCode(alanine, "A").
aminoAcidCode(tyrosine, "Y").
aminoAcidCode(histidine, "H").
aminoAcidCode(glutamine, "Q").
aminoAcidCode(asparagine, "N").
aminoAcidCode(lysine, "K").
aminoAcidCode(asparticAcid, "D").
aminoAcidCode(glutamicAcid, "E").
aminoAcidCode(cysteine, "C").
aminoAcidCode(tryptophan, "W").
aminoAcidCode(arginine, "R").
aminoAcidCode(glycine, "G").

% add amino acids into the kb
aminoAcid(phenylalanine).
aminoAcid(leucine).
aminoAcid(isoleucine).
aminoAcid(methionine).
aminoAcid(valine).
aminoAcid(serine).
aminoAcid(proline).
aminoAcid(threonine).
aminoAcid(alanine).
aminoAcid(tyrosine).
aminoAcid(histidine).
aminoAcid(glutamine).
aminoAcid(asparagine).
aminoAcid(lysine).
aminoAcid(asparticAcid).
aminoAcid(glutamicAcid).
aminoAcid(cysteine).
aminoAcid(tryptophan).
aminoAcid(arginine).
aminoAcid(glycine).

% add protein sequence(s) into the KB (name, sequence)
aminoSequence(spikesite, "KRSFIEDLLFNKV").

aminoAcidPosition(NAME, INDEX, AMINO) :- 
    aminoSequence(NAME, SEQ),
    aminoAcidCode(AMINO, @letteri(SEQ, INDEX)),
    aminoSequenceLength(NAME, LEN),
    INDEX = 1..LEN.

aminoCodeCandidate(NAME, INDEX, AMINO, L1) :-
    aminoAcidPosition(NAME, INDEX, AMINO),
    codon(AMINO, L1).

aminoSequenceLength(NAME, LEN) :- 
    aminoSequence(NAME, STR), 
    LEN = @stringLen(STR).

aminoSequenceCandidateDNASubstring(NAME, INDEX, @concatfun(STR1, STR2)) :-
    INDEX > 1,
    aminoSequenceCandidateDNASubstring(NAME, INDEX-1, STR1),
    aminoCodeCandidate(NAME, INDEX, _, STR2).

aminoSequenceCandidateDNASubstring(NAME, INDEX, STR) :-
    INDEX == 1,
    aminoCodeCandidate(NAME, INDEX, _, STR).

aminoSequenceCandidateDNA(NAME, STR) :-
    aminoSequenceCandidateDNASubstring(NAME, INDEX, STR),
    aminoSequenceLength(NAME, INDEX).

% display answer
seqResult(WHAT):-
    aminoSequenceCandidateDNA(spikesite, WHAT).
#show seqResult/1.