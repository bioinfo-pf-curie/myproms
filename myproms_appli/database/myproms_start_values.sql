/*=============================================================================*/
/*                  myProMS database v3.5.11                                   */
/* MySQL script for required starting values of myProMS v3.6 database          */
/* Requires MySQL 5+                                                           */
/* Patrick Poullet    06/12/2017                                               */
/* Command:                                                                    */
/* >mysql -u <DB_USER> -p -D <DB_NAME> -h <DB_HOST> < myproms_start_values.sql */
/*=============================================================================*/


/*==============================================================*/
/* Loading table USER_PROFILE                                   */
/*==============================================================*/
INSERT INTO USER_PROFILE (ID_PROFILE,NAME) VALUES
(10,'guest'),
(20,'user'),
(21,'power_user'),
(22,'super_user'),
(30,'administrator'),
(31,'power_administrator'),
(32,'super_administrator');


/*==============================================================*/
/* Loading table USER_LIST                                      */
/*==============================================================*/
INSERT INTO USER_LIST (ID_USER,USER_NAME,USER_STATUS,USER_PASSWORD,USER_INFO,USER_PREF) VALUES ('login','No name','bioinfo','Pr4XYhK.C.eos','Default bioinformatician','specApp=1;');


/*==============================================================*/
/* Loading table DATABANK_TYPE                                  */
/*==============================================================*/
INSERT INTO DATABANK_TYPE (ID_DBTYPE,NAME,DES,PARSE_RULES,DEF_IDENT_TYPE,PARSE_STRING,IDENT_MOD_STRING) VALUES
(1,'NCBI - GI','NCBI Non-redundant protein sequences','ID=(gi\\|\\d+\\.?\\d*),:,ORG=\\s+\\[([^\\[\\(]*).*\\]\\Z,:,','GI_ACCESSION','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><B>gi\\|125987826</B><I>\\|sp\\|P15311\\|EZRI_HUMAN</I> ...<BR><B>-Description:</B> >... <B>Ezrin (p81) (Cytovillin) (Villin-2)</B> [...]<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... [<B>Homo sapiens</B> <I>(Human)</I>]',NULL),
(2,'Protein (User-defined) #1','Custom-made list of protein sequences','ID=^(\\S+),:,ORG=(),:,',NULL,'<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><B>All identifiers</B> ...<BR><B>-Description:</B> >... <B>Any string</B><BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>No species</B>',NULL),
(3,'MSDB #1','Non-identical protein sequences','ID=^(\\S+),:,ORG=\\s?- ([^\\(]*),:,',NULL,'<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><B>EZRI_HUMAN</B> ...<BR><B>-Description:</B> >... <B>Ezrin (p81) (Cytovillin) (Villin-2)</B> - ...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... - <B>Homo sapiens</B> <I>(Human)</I>',NULL),
(4,'SWISSPROT/trEMBL #1','Non-identical protein sequences','ID=^(\\S+),:,ORG=\\s?- ([^\\(]*),:,','UNIPROT_ALL','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><B><I>sp\\|</I>P15311\\|EZRI_HUMAN</B> ...<BR><B>-Description:</B> >... <B>Ezrin (p81) (Cytovillin) (Villin-2)</B> - ...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... - <B>Homo sapiens</B> <I>(Human)</I>','.*|->'),
(5,'SWISSPROT/trEMBL #2','Non-identical protein sequences','ID=([A-Z][\\w\\.]+),:,ORG=\\s?- ([^\\(]*),:,','UNIPROT_ACCESSION','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><I>sp\\|</I><B>P15311</B><I>\\|EZRI_HUMAN</I> ...<BR><B>-Description:</B> >... <B>Ezrin (p81) (Cytovillin) (Villin-2)</B> - ...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... - <B>Homo sapiens</B> <I>(Human)</I>',NULL),
(6,'SWISSPROT/trEMBL #3','Non-identical protein sequences','ID=([^\\s\\|]+_[^\\s\\|]+),:,ORG=\\s?- ([^\\(]*),:,','UNIPROT_ID','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><I>sp\\|P15311\\|</I><B>EZRI_HUMAN</B> ...<BR><B>-Description:</B> >... <B>Ezrin (p81) (Cytovillin) (Villin-2)</B> - ...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... - <B>Homo sapiens</B> <I>(Human)</I>',NULL),
(7,'UniProt - ALL','Non-identical protein sequences','ID=^(\\S+),:,DES=^\\S+ (.*)\\sOS=,:,ORG=\\sOS=([^\\(]+?)\\s(\\(|[A-Z]{2}=),:,','UNIPROT_ALL','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><B>sp\\|P15311\\|EZRI_HUMAN</B> ...<BR><B>-Description:</B> >... <B>Ezrin</B> OS=...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... OS=<B>Homo sapiens</B> <I>(Human) GN=...</I>','.*|->'),
(8,'UniProt - ACC','Non-identical protein sequences','ID=([A-Z][\\w\\.-]+),:,DES=^\\S+ (.*)\\sOS=,:,ORG=\\sOS=([^\\(]+?)\\s(\\(|[A-Z]{2}=),:,','UNIPROT_ACCESSION','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><I>sp\\|</I><B>P15311</B><I>\\|EZRI_HUMAN</I> ...<BR><B>-Description:</B> >... <B>Ezrin</B> OS=...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... OS=<B>Homo sapiens</B> <I>(Human) GN=...</I>',NULL),
(9,'UniProt - ID','Non-identical protein sequences','ID=([^\\s\\|]+_[^\\s\\|]+),:,DES=^\\S+ (.*)\\sOS=,:,ORG=\\sOS=([^\\(]+?)\\s(\\(|[A-Z]{2}=),:,','UNIPROT_ID','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><I>sp\\|P15311\\|</I><B>EZRI_HUMAN</B> ...<BR><B>-Description:</B> >... <B>Ezrin</B> OS=...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... OS=<B>Homo sapiens</B> <I>(Human) GN=...</I>',NULL),
(10,'IPI Databank #1','International Protein Index','ID=IPI:([^| .]*),:,DES=^\\S+ (.*),:,ORG=Tax_Id=(\\d+),:,','IPI_ACCESSION','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;>IPI:<B>IPI00843975</B><I>.1\\|SWISS-PROT:P15311\\|...</I><BR><B>-Description:</B> >... <B>Tax_Id=9606 Gene_Symbol=EZR Ezrin</B><BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... Tax_Id=<B>9606</B><I> Gene_Symbol=EZR Ezrin</I>',NULL),
(11,'NCBI - ALL','NCBI Non-redundant protein sequences','ID=^(\\S+),:,ORG=\\s+\\[([^\\[\\(]*).*\\]\\Z,:,','NCBI_ALL','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><B>gi\\|125987826\\|sp\\|P15311\\|EZRI_HUMAN</B> ...<BR><B>-Description:</B> >... <B>Ezrin (p81) (Cytovillin) (Villin-2)</B> [...]<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... [<B>Homo sapiens</B> <I>(Human)</I>]','(gi\\|\\d+)=\$1'),
(12,'UniProt - ACC|ID','Non-identical protein sequences','ID=([A-Z][\\w\\.]+\\|[^\\s\\|]+_[^\\s\\|]+),:,DES=^\\S+ (.*)\\sOS=,:,ORG=\\sOS=([^\\(]+?)\\s(\\(|[A-Z]{2}=),:,','UNIPROT_ALL','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><I>sp\\|</I><B>P15311\\|EZRI_HUMAN</B> ...<BR><B>-Description:</B> >... <B>Ezrin</B> OS=...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... OS=<B>Homo sapiens</B> <I>(Human) GN=...</I>','.*|->'),
(20,'Sigma-Aldrich UPS - ALL','UPS1 & UPS2 Proteomic Standards','ID=^(\\S+),:,ORG=\\s?- ([^\\(]*),:,','UPS_ALL','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><B>P02768ups\\|ALBU_HUMAN_UPS</B> ...<BR><B>-Description:</B> >... <B>Serum albumin (Chain 26-609)</B> - ...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... - <B>Homo sapiens</B> <I>(Human)</I>','.*|->'),
(21,'Sigma-Aldrich UPS - ACC','UPS1 & UPS2 Proteomic Standards','ID=([A-Z][\\w\\.-]+ups),:,ORG=\\s?- ([^\\(]*),:,','UPS_ACCESSION','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;<B>P02768ups</B><I>\\|ALBU_HUMAN_UPS</I> ...<BR><B>-Description:</B> >... <B>Serum albumin (Chain 26-609)</B> - ...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... - <B>Homo sapiens</B> <I>(Human)</I>',NULL),
(22,'Sigma-Aldrich UPS - ID','UPS1 & UPS2 Proteomic Standards','ID=([^\\s\\|]+_[^\\s\\|]+_UPS),:,ORG=\\s?- ([^\\(]*),:,','UPS_ID','<B>-Identifier:</B>&nbsp;&nbsp;&nbsp;&nbsp;><I>P02768ups\\|</I><B>ALBU_HUMAN_UPS</B> ...<BR><B>-Description:</B> >... <B>Serum albumin (Chain 26-609)</B> - ...<BR><B>-Species:</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;>... ... - <B>Homo sapiens</B> <I>(Human)</I>',NULL);



/*==============================================================*/
/* Loading table INSTRUMENT                                     */
/*==============================================================*/
INSERT INTO INSTRUMENT (ID_INSTRUMENT,NAME,DES,COMMENTS,USE_STATUS,UPDATE_DATE,UPDATE_USER,DELTA_PARENT,DELTA_FRAGMENT,TOL_FRAGMENT,NR_LEVEL,RULES) VALUES
(1, 'ESI-FTICR', '', '', 'yes', '2008-02-04 16:30:46', 'myproms', 'ppm', 'Da', 0.8, 0.02, 'b*=1,y°=1,a=2,2+=1,y*=1,y=2,+=1,b=2,b°=1'),
(2, 'MALDI-TOF-TOF', '', '', 'yes', '2008-02-04 17:07:07', 'myproms', 'ppm', 'Da', 0.3, 0.03, 'b*=1,y°=1,a=2,y*=1,ya=2,y=2,+=1,yb=2,b=2,b°=1'),
(3, 'ESI-TRAP', '', '', 'yes', '2008-02-04 16:30:46', 'myproms', 'Da', 'Da', 0.8, 0.02, 'b*=1,y°=1,a=2,2+=1,y*=1,y=2,+=1,b=2,b°=1'),
(4, 'ESI-QUAD-TOF', '', '', 'yes', '2008-02-04 16:30:46', 'myproms', 'Da', 'Da', 0.2, 0.02, 'b*=1,y°=1,a=2,2+=1,y*=1,y=2,+=1,b=2,b°=1'),
(5, 'MALDI-TOF-PSD', '', '', 'yes', '2008-02-04 17:07:07', 'myproms', 'ppm', 'Da', 0.3, 0.03, 'b*=1,y°=1,a=2,y*=1,ya=2,y=2,+=1,yb=2,b=2,b°=1'),
(6, 'ESI-QUAD', '', '', 'yes', '2008-02-04 16:30:46', 'myproms', 'Da', 'Da', 0.6, 0.02, 'b*=1,y°=1,a=2,2+=1,y*=1,y=2,+=1,b=2,b°=1');


/*==============================================================*/
/* Loading table SPECIES                                        */
/*==============================================================*/
INSERT INTO SPECIES (ID_SPECIES,COMMON_NAME,SCIENTIFIC_NAME,TAXONID,IS_REFERENCE) VALUES
(1,'Human','Homo sapiens',9606,1),
(2,'Mouse','Mus musculus',10090,1),
(3,'Fruit fly','Drosophila melanogaster',7227,1),
(4,"Baker's yeast",'Saccharomyces cerevisiae',4932,1),
(5,'E. coli K12','Escherichia coli str. K-12 substr. MG1655',511145,1),
(6,'Caenorhabditis elegans','Caenorhabditis elegans',6239,1);

/*==============================================================*/
/* Loading table IDENTIFIER                                     */
/*==============================================================*/
INSERT INTO IDENTIFIER (ID_IDENTIFIER,ID_SPECIES,NAME,CODE,RESRC_NAME,RESRC_URL,TARGET) VALUES
(1,NULL,'UniProt AC','AC','UniProt Knowledge Base','http://www.uniprot.org/uniprot/XXXX','Protein'),
(2,NULL,'UniProt ID','ID','UniProt Knowledge Base','http://www.uniprot.org/uniprot/XXXX','Protein'),
(3,NULL,'GI Accession','GI','NCBI Protein Entry','http://www.ncbi.nlm.nih.gov/protein/XXXX','Protein'),
(10,NULL,'Gene Name','GN',NULL,NULL,'Gene'),
(20,NULL,'IPI','IPI','Internaltional Protein Index','http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ipi&id=XXXX','Protein'),
(30,NULL,'PIR','PIR','Protein Information Resource','http://pir.georgetown.edu/cgi-bin/nbrfget?uid=XXXX','Protein'),
(40,NULL,'RefSeq','RefSeq','NCBI Reference Sequence','http://www.ncbi.nlm.nih.gov/VVVVXXXX','All'),
(50,NULL,'UniGene','UniGene','UniGene','http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=YYYY&CID=ZZZZ','Gene'),
(60,NULL,'Ensembl','Ensembl','Ensembl Genome Browser','http://www.ensembl.org/id/XXXX','All'),
(70,NULL,'Entrez GeneID','GeneID','Entrez GeneID','http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&term=XXXX','Gene'),
(80,NULL,'KEGG Entry','KEGG','KEGG Pathway Database','http://www.genome.jp/dbget-bin/www_bget?XXXX','Protein'),
(90,NULL,'IntAct','IntAct','IntAct Protein Interaction Database','http://www.ebi.ac.uk/intact/pages/interactions/interactions.xhtml?query=XXXX*','Protein'),
(100,NULL,'STRING','STRING','STRING Protein-Protein Interactions','http://string-db.org/newstring_cgi/show_network_section.pl?identifier=XXXX','Protein'),
(110,4,'SGD','SGD','Saccharomyces Genome Database','http://www.yeastgenome.org/cgi-bin/locus.fpl?dbid=XXXX','Gene'),
(120,3,'FlyBase','FlyBase','FlyBase','http://flybase.org/reports/XXXX.html','Gene'),
(130,2,'MGI','MGI','Mouse Genome Informatics','http://www.informatics.jax.org/marker/XXXX','Gene'),
(140,1,'OMIM','MIM','Online Mendelian Inheritance in Man','http://www.omim.org/entry/XXXX','Gene'),
(150,1,'HPA','HPA','The Human Protein Atlas','http://www.proteinatlas.org/tissue_profile.php?antibody_id=XXXX','Protein'),
(160,6,'WormBase CDS','W_CDS','WormBase','http://www.wormbase.org/db/seq/sequence?name=XXXX;class=CDS','CDS'),
(161,6,'WormBase Protein','W_P','WormBase','http://www.wormbase.org/db/seq/protein?name=WP:XXXX;class=Protein','Protein'),
(162,6,'WormBase Gene','W_G','WormBase','http://www.wormbase.org/db/gene/gene?name=XXXX;class=Gene','Gene');

/*==============================================================*/
/* Loading table MODIFICATION                                   */
/*==============================================================*/
INSERT INTO MODIFICATION (ID_MODIFICATION,UNIMOD_ACC,PSI_MS_NAME,INTERIM_NAME,DES,SYNONYMES,COMPOSITION,MONO_MASS,AVGE_MASS,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS) VALUES
(1,35,'Oxidation','Hydroxylation','Oxidation or Hydroxylation','','O',15.9949,15.9994,'C,D,F,G;*,H,K,M,N,P,R,W,Y','O','FF0000',0,0,1),
(2,1,'Acetyl','Acetyl','Acetylation','','H(2) C(2) O',42.0106,42.0367,'S,T,=,K,Y,-,H,C','A','00CC00',0,0,1),
(3,4,'Carbamidomethyl','Carbamidomethyl','Iodoacetamide derivative','##Carboxyamidomethylation##Iodoacetamide##Iodoacetamide##','H(3) C(2) N O',57.0215,57.0513,'S,T,=,K,Y,E,H,D,C','C','FF9900',0,0,1),
(4,58,'Propionyl','Propionyl_light','Propionate labeling reagent light form (N-term & K)','','H(4) C(3) O',56.0262,56.0633,'-,S,T,=,K','P','00CC00',0,0,1),
(5,34,'Methyl','Methyl','Methylation','##methyl ester##','H(2) C',14.0156,14.0266,'S,T,N,=,K,*,E,H,-,Q,C,D,I,R,L','M','0000FF',0,0,1),
(6,37,'Trimethyl','tri-Methylation','tri-Methylation','','H(6) C(3)',42.047,42.0797,'A;-,K,R','T','0000FF',0,0,1),
(7,121,'GG','GlyGly','ubiquitinylation residue','##glycineglycine##','H(6) C(4) N(2) O(2)',114.043,114.103,'S,T,C,K','U','FF0000',0,0,1),
(8,36,'Dimethyl','di-Methylation','di-Methylation','','H(4) C(2)',28.0313,28.0532,'K,N,P;-,R,=','D','0000FF',0,0,1),
(9,21,'Phospho','Phospho','Phosphorylation','','H O(3) P',79.9663,79.9799,'C,D,H,K,R,S,T,Y','P','FF9900',0,0,1);


/*==============================================================*/
/* Loading table QUANTIFICATION_METHOD                          */
/*==============================================================*/
INSERT INTO QUANTIFICATION_METHOD (ID_QUANTIFICATION_METHOD,NAME,CODE,DES) VALUES
(1,'Ext. ion chrom.','XIC','Extracted ion chromatogram'),
(2,'Prot. ratio pep.','PROT_RATIO_PEP','Protein ratio based on peptide ratio'),
(3,'SILAC','SILAC','Stable isotope labeling with amino acids in cell culture'),
(4,'iTRAQ','ITRAQ','Isobaric tag for relative and absolute quantitationStable isotope labeling'),
(5,'emPAI','EMPAI','Exponentially Modified Protein Abundance Index'),
(6,'SIN','SIN','Normalized Spectral Index'),
(7,'TnPQ','TNPQ','Method based on T3PQ <-> Top 3 Protein Quantification'),
(8,'MaxQuant','MQ','MaxQuant related Quantitation (iBAQ, Intensity, LFQ and MS/MS count)'),
(9,'DIA','DIA','Data independant acquisition'),
(10,'SSP Analysis','SSPA','State-Specific Protein Analysis'),
(11,'TMT','TMT','Tandem mass tag');

/*==============================================================*/
/* Loading table QUANTIFICATION_PARAMETER                       */
/*==============================================================*/
INSERT INTO QUANTIFICATION_PARAMETER (ID_QUANTIF_PARAMETER,ID_QUANTIFICATION_METHOD,NAME,CODE,DES) VALUES
(1,1,'XIC intensity','XIC_INTENSITY','Intensity of extracted ion chromatogram'),
(2,1,'XIC area','XIC_AREA','Area of extracted ion chromatogram'),
(67,1,'Retention-time begin','RT_BEGIN','Starting chromatogram point to compute XIC area'),
(68,1,'Retention-time end','RT_END','Ending chromatogram point to compute XIC area'),
(69,1,'Retention time apex','RT_APEX','Retention time at pic apex'),

(4,3,'Isotope intensity','ISO_INTENSITY','Intensity of isotope ion chromatogram'),
(5,3,'Isotope area','ISO_AREA','Area of isotope ion chromatogram'),
(71,3,'Retention time apex','RT_APEX','Retention time at pic apex'),
(72,3,'Retention-time begin','RT_BEGIN','Starting chromatogram point to compute XIC area'),
(73,3,'Retention-time end','RT_END','Ending chromatogram point to compute XIC area'),

(7,4,'Reporter intensity','REP_INTENSITY','Intensity of reporter ion chromatogram'),
(8,4,'Reporter area','REP_AREA','Area of reporter ion chromatogram'),
(9,4,'Reporter mass','REP_MASS','Observed mass of reporter ion'),

(11,2,'p-value','PVAL','p-value for protein ratio'),
(12,2,'Adj. p-value','PVAL_ADJ','Adjusted p-value for protein ratio'),
(13,2,'Ratio','RATIO',NULL),
(14,2,'Geometric std. dev.','SD_GEO',NULL),
(15,2,'Low. bound for conf. Int','CONF_LOW',NULL),
(16,2,'Up. bound for conf. Int','CONF_UP',NULL),
(17,2,'# pep. used','NUM_PEP_USED','Number of peptides used for protein ratio calculation'),
(18,2,'Normality p-value','NORM_PVAL','p-value for normality test of peptide ratios distribution'),
(19,2,'Adj. Normality p-value','NORM_PVAL_ADJ','Adjusted p-value for normality test of peptide ratios distribution'),
(20,2,'T-test power','TTEST_POWER','Power of T-test used when peptide ratios distribution is normal'),
/*(21,2,'Sign. change','SIG_CHANGE','Change in protein abundance is significant'),*/
(22,2,'# pep. total','NUM_PEP_TOTAL','Number of peptides provided for protein ratio calculation'),
(23,2,'Mean repl.','MEAN_STATE','Average of replicates for selected state'),
(24,2,'Variation coeff.','COEF_VAR_STATE','Variation coefficient between replicates for selected state'),
(25,2,'Var. coeff. exclusion','COEF_VAR_EXCL','Exclusion based on coefficient of variation'),
(26,2,'Ratio exclusion','RATIO_EXCL','Exclusion based on ratio p-value'),
(27,2,'Excluded','EXCL','Excluded from statistics'),
(28,2,'# dist. pep. used','DIST_PEP_USED','Number of distinct peptides used for protein ratio calculation'),
(81,2,'Peptides','PEPTIDES','Number of peptides sequences'),
(82,2,'Razor + unique','RAZ_UNI_PEP','Number of razor + unique peptides'),
(83,2,'Unique','UNIQUE_PEP','Number of unique peptides'),
(84,2,'Ratio','RATIO_RAW','Non-normalized ratio'),
(85,2,'Ratio variability','RATIO_VAR','Coefficient of variability over all redundant quantifiable peptides (%)'),

(31,5,'emPAI','EMPAI','emPAI as computed by Mascot web server'),
(32,5,'emPAI mol. percent.','EMPAI_MOL','emPAI molar fraction percentage (emPAI divided by tot. sum of emPAI)'),
(33,5,'emPAI Mr percent.','EMPAI_MR','emPAI weight fraction percentage (emPAI divided by tot. sum of emPAIxMr)'),

(41,6,'SI of a peptide','SIN_SI','Sum for a peptide sequence of all the intensities of the queries that match this pepseq'),
(42,6,'SIN','SIN_SIN','Spectral Index Normalized (sum of SI of all peptides of a protein divided by SIGI (global intensity)'),

(51,7,'Mean TnPQ','TNPQ_MEAN','Mean of the 3 best peptides intensity identified in the protein'),
(52,7,'p-value','PVAL','p-value for protein ratio'),
(53,7,'Adj. p-value','PVAL_ADJ','Adjusted p-value for protein ratio'),
(54,7,'Ratio','RATIO','ratio of mean intensities between two conditions'),
(55,7,'Low. bound for conf. Int','CONF_LOW',NULL),
(56,7,'Up. bound for conf. Int','CONF_UP',NULL),
(57,7,'Normality p-val.','NORM_PVAL','p-value threshold for normality test of peptide ratios distribution'),
(58,7,'Adj. normality p-val.','NORM_PVAL_ADJ','Adjusted p-value threshold for normality test of peptide ratios distribution'),
(59,7,'T-test power','TTEST_POWER','Power of T-test used when peptide ratios distribution is normal'),
(60,7,'# pep. used','NUM_PEP_USED','Number of peptides used for protein ratio calculation'),
(61,7,'# pep. total','NUM_PEP_TOTAL','Number of peptides provided for protein ratio calculation'),
(62,7,'# dist. pep. used','DIST_PEP_USED','Number of distinct peptides used for protein ratio calculation'),

(91,8,'Intensity','MQ_INT','Intensity column of MaxQuant output program'),
(92,8,'iBAQ','MQ_IBAQ','intensity-based absolute quantification'),
(93,8,'LFQ','MQ_LFQ','maxLFQ or label-free quantification'),
(94,8,'MS/MS count','MQ_SC','Spectral count'),
(95,8,'Peptides','PEPTIDES','Number of peptides sequences'),
(96,8,'Razor + unique','RAZ_UNI_PEP','Number of razor + unique peptides'),
(97,8,'Unique','UNIQUE_PEP','Number of unique peptides'),

(101,10,'State mean','MEAN','Mean of peptide count of all state replicates'),
(102,10,'Effect','EFFECT','Log of difference between sets of states'),
(103,10,'Best delta','DELTA_PC','Largest mean delta between states normalized by best mean (in %)'),
(104,10,'p-value','PVAL','p-value univariate test'),
(105,10,'Adj. p-value','PVAL_ADJ','Adjusted p-value of univariate test'),

(111,11,'Reporter intensity','REP_INTENSITY','Intensity of reporter ion chromatogram'),
(112,11,'Reporter area','REP_AREA','Area of reporter ion chromatogram'),
(113,11,'Reporter mass','REP_MASS','Observed mass of reporter ion');


/*==============================================================*/
/* Loading table REFERENCE_RT                                   */
/*==============================================================*/
INSERT INTO REFERENCE_RT (ID_REFERENCE_RT,NAME,DES,COLUMN_TYPE,PROT_IDENTIFIER,PROT_SEQ,NUM_PEP,DATA) VALUES
(
	1,
	'iRT-C18',
	'Biognosis indexed Retention Time (iRT) on C18 column',
	'C18',
	'Biognosys|iRT Kit Fusion',
	'AGGSSEPVTGLADKVEATFGVDESANKYILAGVESNKDAVTPADFSEWSKFLLQFGAQGSPLFKLGGNETQVRTPVISGGPYYERTPVITGAPYYERGDLDAASYYAPVRTGFIIDPGGVIRGTFIIDPAAIVR',
	11,
'<PEPTIDE_DATA version="1.0.1">
<PEPTIDE pepId="1" sequence="LGGNETQVR" monoMass="487.257" charge="2" iRT="-24.92"/>
<PEPTIDE pepId="2" sequence="AGGSSEPVTGLADK" monoMass="644.823" charge="2" iRT="0.00"/>
<PEPTIDE pepId="3" sequence="VEATFGVDESANK" monoMass="683.828" charge="2" iRT="12.39"/>
<PEPTIDE pepId="4" sequence="YILAGVESNK" monoMass="547.298" charge="2" iRT="19.79"/>
<PEPTIDE pepId="5" sequence="TPVISGGPYYER" monoMass="669.838" charge="2" iRT="28.71"/>
<PEPTIDE pepId="6" sequence="TPVITGAPYYER" monoMass="683.854" charge="2" iRT="33.38"/>
<PEPTIDE pepId="7" sequence="GDLDAASYYAPVR" monoMass="699.339" charge="2" iRT="42.26" excluded="1"/>
<PEPTIDE pepId="8" sequence="DAVTPADFSEWSK" monoMass="726.836" charge="2" iRT="54.62"/>
<PEPTIDE pepId="9" sequence="TGFIIDPGGVIR" monoMass="622.854" charge="2" iRT="70.52"/>
<PEPTIDE pepId="10" sequence="GTFIIDPAAIVR" monoMass="636.869" charge="2" iRT="87.23"/>
<PEPTIDE pepId="11" sequence="FLLQFGAQGSPLFK" monoMass="776.93" charge="2" iRT="100.00"/>
</PEPTIDE_DATA>'
),
(
	3,
	'iRT-C18 v2',
	'Biognosis indexed Retention Time (iRT) on C18 column. v2',
	'C18',
	'Biognosys|iRT-Kit_WR_fusion',
	'LGGNEQVTRYILAGVENSKGTFIIDPGGVIRGTFIIDPAAVIRGAGSSEPVTGLDAKTPVISGGPYEYRVEATFGVDESNAKTPVITGAPYEYRDGLDAASYYAPVRADVTPADFSEWSKLFLQFGAQGSPFLK',
	11,
'<PEPTIDE_DATA version="1.0.1">
<PEPTIDE pepId="1" sequence="LGGNEQVTR" monoMass="487.257" charge="2" iRT="-24.92"/>
<PEPTIDE pepId="2" sequence="GAGSSEPVTGLDAK" monoMass="644.823" charge="2" iRT="0.00"/>
<PEPTIDE pepId="3" sequence="VEATFGVDESNAK" monoMass="683.828" charge="2" iRT="12.39"/>
<PEPTIDE pepId="4" sequence="YILAGVENSK" monoMass="547.298" charge="2" iRT="19.79"/>
<PEPTIDE pepId="5" sequence="TPVISGGPYEYR" monoMass="669.838" charge="2" iRT="28.71"/>
<PEPTIDE pepId="6" sequence="TPVITGAPYEYR" monoMass="683.854" charge="2" iRT="33.38"/>
<PEPTIDE pepId="7" sequence="DGLDAASYYAPVR" monoMass="699.339" charge="2" iRT="42.26" excluded="1"/>
<PEPTIDE pepId="8" sequence="ADVTPADFSEWSK" monoMass="726.836" charge="2" iRT="54.62"/>
<PEPTIDE pepId="9" sequence="GTFIIDPGGVIR" monoMass="622.854" charge="2" iRT="70.52"/>
<PEPTIDE pepId="10" sequence="GTFIIDPAAVIR" monoMass="636.869" charge="2" iRT="87.23"/>
<PEPTIDE pepId="11" sequence="LFLQFGAQGSPFLK" monoMass="776.93" charge="2" iRT="100.00"/>
</PEPTIDE_DATA>'
);


/*==============================================================*/
/* Commit all data                                              */
/*==============================================================*/
COMMIT;
