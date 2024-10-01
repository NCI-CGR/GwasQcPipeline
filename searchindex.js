Search.setIndex({docnames:["api/cgr_gwas_qc","api/config","api/parsers","dev_docs/documentation","dev_docs/setup","dev_docs/source_code","dev_docs/testing","getting_started/configuration","getting_started/installation","getting_started/running_pipeline","index","reference/file_types","sub_workflows/contamination","sub_workflows/delivery","sub_workflows/entry_points","sub_workflows/intensity_check","sub_workflows/sample_qc","sub_workflows/subject_qc"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":3,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":2,"sphinx.domains.rst":2,"sphinx.domains.std":1,"sphinx.ext.todo":2,"sphinx.ext.viewcode":1,sphinx:56},filenames:["api/cgr_gwas_qc.rst","api/config.rst","api/parsers.rst","dev_docs/documentation.rst","dev_docs/setup.rst","dev_docs/source_code.rst","dev_docs/testing.rst","getting_started/configuration.rst","getting_started/installation.rst","getting_started/running_pipeline.rst","index.rst","reference/file_types.rst","sub_workflows/contamination.rst","sub_workflows/delivery.rst","sub_workflows/entry_points.rst","sub_workflows/intensity_check.rst","sub_workflows/sample_qc.rst","sub_workflows/subject_qc.rst"],objects:{"cgr-config":{"--bpm-file":[9,5,1,"cmdoption-cgr-config-bpm-file"],"--cgems":[9,5,1,"cmdoption-cgr-config-cgems"],"--cgems-dev":[9,5,1,"cmdoption-cgr-config-cgems-dev"],"--genome-build":[9,5,1,"cmdoption-cgr-config-genome-build"],"--include-unused-settings":[9,5,1,"cmdoption-cgr-config-u"],"--project-name":[9,5,1,"cmdoption-cgr-config-project-name"],"--sample-sheet":[9,5,1,"cmdoption-cgr-config-s"],"--slurm-partition":[9,5,1,"cmdoption-cgr-config-slurm-partition"],"-s":[9,5,1,"cmdoption-cgr-config-s"],"-u":[9,5,1,"cmdoption-cgr-config-u"],"-y":[9,5,1,"cmdoption-cgr-config-y"]},"cgr-pre-flight":{"--cluster-group-size":[9,5,1,"cmdoption-cgr-pre-flight-cluster-group-size"],"--config-file":[9,5,1,"cmdoption-cgr-pre-flight-config-file"],"--no-reference-files-check":[9,5,1,"cmdoption-cgr-pre-flight-no-reference-files-check"],"--no-update-config":[9,5,1,"cmdoption-cgr-pre-flight-no-update-config"],"--no-user-files-check":[9,5,1,"cmdoption-cgr-pre-flight-no-user-files-check"],"--threads":[9,5,1,"cmdoption-cgr-pre-flight-j"],"-j":[9,5,1,"cmdoption-cgr-pre-flight-j"]},"cgr-submit":{"--biowulf":[9,5,1,"cmdoption-cgr-submit-biowulf"],"--ccad2":[9,5,1,"cmdoption-cgr-submit-ccad2"],"--cgems":[9,5,1,"cmdoption-cgr-submit-cgems"],"--cluster-profile":[9,5,1,"cmdoption-cgr-submit-cluster-profile"],"--dry-run":[9,5,1,"cmdoption-cgr-submit-dry-run"],"--local-mem-mb":[9,5,1,"cmdoption-cgr-submit-local-mem-mb"],"--local-tasks":[9,5,1,"cmdoption-cgr-submit-local-tasks"],"--no-biowulf":[9,5,1,"cmdoption-cgr-submit-biowulf"],"--no-ccad2":[9,5,1,"cmdoption-cgr-submit-ccad2"],"--no-cgems":[9,5,1,"cmdoption-cgr-submit-cgems"],"--no-dry-run":[9,5,1,"cmdoption-cgr-submit-dry-run"],"--no-notemp":[9,5,1,"cmdoption-cgr-submit-notemp"],"--no-slurm":[9,5,1,"cmdoption-cgr-submit-slurm"],"--notemp":[9,5,1,"cmdoption-cgr-submit-notemp"],"--queue":[9,5,1,"cmdoption-cgr-submit-queue"],"--slurm":[9,5,1,"cmdoption-cgr-submit-slurm"],"--submission-cmd":[9,5,1,"cmdoption-cgr-submit-submission-cmd"],"--subworkflow":[9,5,1,"cmdoption-cgr-submit-subworkflow"],"--time-hr":[9,5,1,"cmdoption-cgr-submit-time-hr"]},"cgr_gwas_qc.parsers":{bim:[2,0,0,"-"],bpm:[2,0,0,"-"],eigensoft:[2,0,0,"-"],graf:[2,0,0,"-"],king:[2,0,0,"-"],plink:[2,0,0,"-"]},"cgr_gwas_qc.parsers.bim":{BimFile:[2,1,1,""],BimRecord:[2,1,1,""],open:[2,3,1,""]},"cgr_gwas_qc.parsers.bim.BimRecord":{get_record_problems:[2,2,1,""]},"cgr_gwas_qc.parsers.bpm":{BpmFile:[2,1,1,""],BpmRecord:[2,1,1,""],open:[2,3,1,""]},"cgr_gwas_qc.parsers.bpm.BpmFile":{write:[2,2,1,""]},"cgr_gwas_qc.parsers.eigensoft":{Eigenvec:[2,1,1,""]},"cgr_gwas_qc.parsers.eigensoft.Eigenvec":{components:[2,4,1,""],filename:[2,4,1,""],values:[2,4,1,""]},"cgr_gwas_qc.parsers.graf":{read_relatedness:[2,3,1,""]},"cgr_gwas_qc.parsers.illumina":{IlluminaBeadArrayFiles:[2,0,0,"-"],adpc:[2,0,0,"-"]},"cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles":{BeadPoolManifest:[2,1,1,""],GenotypeCalls:[2,1,1,""],LocusEntry:[2,1,1,""],complement:[2,3,1,""],read_byte:[2,3,1,""],read_char:[2,3,1,""],read_float:[2,3,1,""],read_int:[2,3,1,""],read_scanner_data:[2,3,1,""],read_string:[2,3,1,""],read_ushort:[2,3,1,""]},"cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles.BeadPoolManifest":{addresses:[2,4,1,""],chroms:[2,4,1,""],control_config:[2,4,1,""],manifest_name:[2,4,1,""],normalization_lookups:[2,4,1,""],num_loci:[2,4,1,""],ref_strands:[2,4,1,""],snps:[2,4,1,""],source_strands:[2,4,1,""]},"cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles.GenotypeCalls":{get_autocall_date:[2,2,1,""],get_autocall_version:[2,2,1,""],get_ballele_freqs:[2,2,1,""],get_base_calls:[2,2,1,""],get_base_calls_forward_strand:[2,2,1,""],get_base_calls_generic:[2,2,1,""],get_base_calls_plus_strand:[2,2,1,""],get_call_rate:[2,2,1,""],get_cluster_file:[2,2,1,""],get_control_x_intensities:[2,2,1,""],get_control_y_intensities:[2,2,1,""],get_gc10:[2,2,1,""],get_gc50:[2,2,1,""],get_gender:[2,2,1,""],get_genotype_scores:[2,2,1,""],get_genotypes:[2,2,1,""],get_imaging_date:[2,2,1,""],get_logr_dev:[2,2,1,""],get_logr_ratios:[2,2,1,""],get_normalization_transforms:[2,2,1,""],get_normalized_intensities:[2,2,1,""],get_num_calls:[2,2,1,""],get_num_intensity_only:[2,2,1,""],get_num_no_calls:[2,2,1,""],get_num_snps:[2,2,1,""],get_percentiles_x:[2,2,1,""],get_percentiles_y:[2,2,1,""],get_ploidy:[2,2,1,""],get_ploidy_type:[2,2,1,""],get_raw_x_intensities:[2,2,1,""],get_raw_y_intensities:[2,2,1,""],get_sample_name:[2,2,1,""],get_sample_plate:[2,2,1,""],get_sample_well:[2,2,1,""],get_scanner_data:[2,2,1,""],get_slide_identifier:[2,2,1,""],get_snp_manifest:[2,2,1,""],is_write_complete:[2,2,1,""],supported_versions:[2,4,1,""]},"cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles.LocusEntry":{address_a:[2,4,1,""],address_b:[2,4,1,""],assay_type:[2,4,1,""],chrom:[2,4,1,""],map_info:[2,4,1,""],name:[2,4,1,""],ref_strand:[2,4,1,""],snp:[2,4,1,""],source_strand:[2,4,1,""]},"cgr_gwas_qc.parsers.illumina.adpc":{AdpcBase:[2,1,1,""],AdpcReader:[2,1,1,""],AdpcWriter:[2,1,1,""]},"cgr_gwas_qc.parsers.king":{read_related:[2,3,1,""]},"cgr_gwas_qc.parsers.plink":{read_genome:[2,3,1,""],read_het:[2,3,1,""],read_hwe:[2,3,1,""],read_imiss:[2,3,1,""],read_lmiss:[2,3,1,""],read_sexcheck:[2,3,1,""]},"cgr_gwas_qc.testing.data":{FakeData:[6,1,1,""],RealData:[6,1,1,""]},"cgr_gwas_qc.workflow.scripts":{agg_population_qc_tables:[17,0,0,"-"],sample_qc_table:[16,0,0,"-"],snp_qc_table:[16,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","function","Python function"],"4":["py","attribute","Python attribute"],"5":["std","cmdoption","program option"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:function","4":"py:attribute","5":"std:cmdoption"},terms:{"00003":2,"001":9,"001_1_0000000":7,"001_1_analysismanifest_0000000":7,"0_20011747_a1":[7,11],"0x01":11,"0x0f":11,"0x1b":11,"0x6b":11,"0x6c":11,"0xdc":11,"0xe7":11,"100":[2,9],"1000":[7,9,12,16],"1000genom":11,"10th":[2,11],"11100111":11,"120":9,"12th":11,"15_grch38_full_analysis_set":7,"1st":[2,7],"2010":2,"2012":11,"20130502":7,"2014":2,"2015":2,"2017":[2,11],"20240227130627":7,"20241001152427":7,"20g":8,"231":9,"24v1":[7,11],"2867":2,"2873":2,"2nd":[2,7],"3rd":[2,5,7,11],"4th":11,"50k":9,"50th":2,"5th":[2,11],"6000":7,"6336":7,"6th":11,"700078":7,"7th":11,"8000":9,"8gb":9,"8th":11,"9297":2,"95th":2,"9th":11,"boolean":[7,16],"byte":[2,11],"case":[2,6,7,9,11,16,17],"catch":4,"char":2,"class":[1,2,5,6],"default":[1,4,5,7,8,9,12,13,14,15,16,17],"export":4,"final":[9,11,16,17],"float":[2,16,17],"function":[2,5,6,11,16],"import":[1,5,6],"int":[2,6,17],"long":7,"new":[3,4,5,9,16],"null":[2,4],"return":[1,2,5,6],"sch\u00e4ffer":2,"static":[2,3],"switch":7,"throw":9,"true":[2,3,6,7,16],"try":5,"while":[2,4,5,6,9,10],Age:11,And:[5,8],For:[2,3,4,5,6,7,9,10,11,16],IBS:[2,11],IDs:[2,11,17],Ids:7,NOT:2,Not:9,One:[2,5,11],PCs:11,The:[1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17],Then:[3,6,9,17],There:[4,5,6,7,9,11,14],These:[4,5,7,8,9,11],Use:[1,5,9],Useful:2,Using:5,With:9,__init__:[5,6],__main__:5,_build:[3,4],_grn:7,_ld:17,_pb:11,_red:7,_snps_autosom:17,abf:12,about:[2,5,7,9,11],abov:[3,7,11],absolut:[1,4,5],accept:7,access:[4,5,8],accompani:11,across:16,actgid:11,action:3,activ:[1,4,5,8],actual:7,add:[2,3,4,5,6,9,16,17],add_user_fil:6,added:[6,9],adding:[2,5,9,10],addit:[2,4,7,9,17],address:2,address_a:2,address_b:2,addressa:2,addressa_id:11,addressb:2,addressb_id:11,adjust:[3,5,9],adpc:[2,5],adpcbas:2,adpcread:2,adpcwrit:2,advanc:3,advantag:[9,10],aff:[2,11],afr:16,afr_af:7,african:16,after:[7,8,11,16],ag_match:2,ag_miss:2,again:[5,9,17],agg_contamin:7,agg_population_qc_t:17,aggreg:[7,12,14,15,16,17],agmr:2,ahead:4,ajhg:2,alk:2,all:[1,2,4,5,6,7,9,10,11,16,17],all_sample_qc:13,allel:[2,7,11,12],allele_1:2,allele_2:2,allelea_probeseq:11,alleleb_probeseq:11,alloc:7,allow:[2,3,5,7,9],almost:[5,9],along:9,alpha:11,alphanumer:2,alreadi:[3,6,7],also:[3,4,5,6,7,9,11,13,16,17],alwai:[9,11],amaz:10,ambiguous_allel:2,among:6,amount:[5,6,9,11],amr_af:7,anaconda:[4,8],analysi:7,analysismanifest:13,analytic_exclus:16,analytic_exclusion_reason:16,ancestr:[10,11,16,17],ancestri:[2,11],ancestry_s1:11,ancestry_s2:11,ancestry_s3:11,ani:[4,6,7,8,9,11,16],annoi:9,annot:2,anoth:[3,7,11],anyth:7,api:[3,5,10],appear:4,applic:5,approach:[6,16],arbitrari:[2,17],architectur:5,arg:9,argument:[2,9],around:9,arrai:[2,7,11],array_snp_id:16,artifact:4,asi:11,asian:16,asn:16,aspect:5,aspx:2,assai:2,assay_typ:2,assert:6,asset:4,assign:[2,16,17],associ:[2,11],assum:[1,2,5,7,11],assumpt:5,attempt:7,attribut:2,auto:[3,4],autocal:[2,11],autoconvert:[2,11],automat:[2,3,4,6,7,9],automodul:3,autosom:[2,7],autosomal_het_threshold:[7,17],autosomal_heterozygosity_plot:17,avail:[6,7,9],averag:11,backend:4,bad_posit:2,baf:11,base:[2,4,8,11,16,17],basecal:2,bash:[4,8],bashrc:4,basic:[5,7],basic_stat:2,bcf:[7,12,14,15],bead:[2,11],beadchip:11,beadpoolmanifest:2,beadsetid:11,becaus:[2,4,5,6,9],bed:[6,7,13,14,16,17],befor:[1,5,7,9,16],behavior:[7,16],below:[4,7,9,16],benefit:4,besid:5,between:7,beyond:11,bgz:7,bgzip:[5,7],biallel:11,bim:[2,5,6,7,13,14,16,17],bimfil:2,bimrecord:2,bin:[2,3,4,8],binari:[2,5,11],binomi:[2,11],bioconda:4,bioinformat:[2,6],biowulf:[5,9],bit:[5,11],black:4,blob:[12,13,14,15,16,17],block:[2,11],bold:7,bool:[2,6],bot:11,both:[4,5,6,11],box:[4,7],bpm2abf:7,bpm:[2,5,7,9,12,16],bpm_file:9,bpmfile:2,bpmrecord:2,briefli:[5,11],broadinstitut:2,broken:7,build:[4,5,7,9,11,16,17],built:[3,4,10],bulletin:11,cach:6,calcul:[2,5,7,15,17],call:[2,4,5,6,7,17],call_rat:16,call_rate_1:16,call_rate_2:[16,17],call_rate_initi:16,callabl:[1,5],can:[2,4,5,6,7,9,10,11,14],cannot:11,capabl:6,case_control:[7,9,16,17],case_control_column:[7,9],case_control_dtyp:[16,17],case_control_gwa:7,categor:16,categori:[2,16],ccad2:9,ccad:[5,6,7,9],cell:2,centimorgan:11,central:[5,7],ceu:11,cfg:[1,5],cgem:[5,6,7,9],cgems_copi:5,cgems_jobscript:5,cgems_statu:5,cgems_submit:5,cgf:[5,6,7],cgr:[3,4,5,6,8,10,12,13,14,15,16,17],cgr_gwas_qc:[1,2,4,5,6,8,12,13,14,15,16,17],cgr_refer:5,cgr_sample_sheet:[1,5,7,9],cgroup1:[1,5],cgroup2:[1,5],cgroup:9,challeng:6,chang:[4,5,7,16],channel:4,charact:[2,11],chb:11,chd:11,check:[2,4,5,7,10,16,17],check_sex:2,check_write_complet:2,checker:4,chen:2,chose:5,chr:[2,7,11],chrom:2,chromosom:[2,7,11,16],chrx_inbreed:16,chry:11,classmethod:[1,5],clean:5,cleaner:4,clear:11,cli:[6,9],click:3,clinvar_acmg:11,clone:[4,6],close:[2,3,6],clst:11,cluster:[1,2,7,10,11],cluster_group:[1,5,9],cluster_group_s:9,cluster_profil:[5,6,9],cmd:9,code2genotyp:2,code:[1,2,3,4,10,11],codespel:4,coeffici:[2,7,11,16,17],cog:[2,11],column:[1,2,5,7,9,11],colun:7,com:[2,4,6,8,9,11,12,13,14,15,16,17],combin:[1,2,5],combinatori:[1,5],command:9,comment:16,commit:3,common:[2,5,9,16],common_url:3,commonli:5,compar:[5,11],comparison:16,compat:[1,5],complement:2,complet:[2,9,10,16],complic:[5,11],compon:[2,16,17],componenet:2,comput:[6,11,17],concaten:[9,16,17],concept:6,concord:[2,7,17],concordance_t:7,conda:[1,5,8,9],conda_dir:[1,5],conf:3,confid:11,config:[1,4,5,6,10,11,12,13,14,15,16,17],config_fil:9,configmgr:[1,5],configur:[0,3,10],confpath:3,conftest:6,consid:[2,7,11],consist:5,constant:5,contain:[1,3,5,6,7,11,17],contam:7,contam_popul:[7,12],contam_threshold:[7,16],contamin:[5,7,9,10,11,15],contamination_r:16,content:16,context:2,contig:11,continent:16,continu:9,control:[2,5,7,9,11,16,17],control_config:2,control_hwp_threshold:[7,17],control_statu:[7,11],controls_unrelated_maf:17,conveni:[2,9],convent:6,convert:[5,7,14],coordin:11,copi:[6,7,9,11,13],core:[2,5,9],correct:[5,6,7,13],correctli:[3,7],correspond:[5,7,11],could:6,count:[2,5,11],cpu:9,creat:[1,3,5,6,7,8,14,16],creation:9,criteria:[10,16],csv:[1,5,7,9,12,13,15,16,17],ctrl:2,curl:4,current:[1,2,5,6,7,8,9],custom:[5,9,10],cutoff:7,dai:11,dali:2,data:[1,2,4,5,7,9,11,16],data_dictionari:1,databas:2,datafram:[1,2,5],dataset:[6,11],date:[2,7],dceg:[6,7],dd_dir:1,deactiv:8,decemb:2,decid:[5,7,9],decor:6,def:6,defin:[5,6,9],definit:5,defq:7,degre:[2,7],delet:9,deliver:7,deliveri:[5,7,9,10],depend:7,deploi:10,describ:[5,7,10,11,16],descript:[2,5,8,11,16,17],design:[2,5,9,10,11,13],detail:[2,3,5,7,8,9,11,12,13,14,15,16,17],detect:4,determin:[2,16],dev:[3,4,7,9],dev_doc:3,develop:[3,6,10],deviat:[2,11],dictioanri:1,dictionari:5,did:[9,16],differ:[2,5,6,7,9,11],diploid:11,dir:[6,7],direct:5,directli:[2,3,5,6,9],directori:[1,5,6,7,9],disabl:3,discord:7,discuss:5,dist:4,distanc:[2,11],distinguish:11,distribut:[6,7],divid:[2,10,11],divis:11,dna:11,doc:[3,4,5,9],document:[4,9,11],docx:[1,5,13],docx_templ:[1,5],doe:[6,9,11],don:[2,9],done:7,doubl:2,down:[4,6,7],download:[5,7,8,11],drag:4,driven:6,drop:[4,7],dry:9,dst:[2,11],dtype:[2,16,17],dummi:11,dup:7,dup_concordance_cutoff:[7,16,17],duplic:2,dure:[2,4,7,9,16,17],dynam:5,e0179106:2,e_het:2,e_hom:[2,17],each:[2,4,5,7,9,10,11,12,15,16,17],eas_af:7,easi:[3,4,5,10],easier:[2,5,16],easiest:7,easili:9,east:16,echo:4,edit:[4,7,9],edu:2,egt:7,eigensoft:[0,5,17],eigenvalu:2,eigenvec:[2,17],eighth:11,either:[2,11,13,16],elet:2,els:5,emb:[3,5],embed:4,enabl:4,enchar:2,encod:[2,7,11],encoded_chrom:2,end:7,endchar:2,enhanc:5,ensur:[1,4,6,9],entir:[5,7,9],entri:[2,4,5,7,10],entry_point:[5,9,14],env:[1,3,4,5,8],environ:[3,5,6,8,9,10],environment:6,equilibrium:[2,11,17],error:[4,5,9,10],especi:9,essenti:5,estim:[2,7,12,16,17],estmat:2,etc:11,eur:16,eur_af:7,european:16,even:4,everi:4,everyth:[4,7],exact:[2,11],exactli:7,examin:10,exampl:[1,2,3,4,5,6,9,10,11],exce:16,except:[5,6],exclud:[7,9,11,16],exclus:[5,7],exclusivemaximum:7,exclusiveminimum:7,execut:[7,9],executablepath:3,exist:[6,9],existing_fil:6,exp_clust:11,expampl:7,expand:[1,5],expans:[1,5],expect:[2,5,7,9,11,16,17],expected_sex:[7,9,11,16],expected_sex_column:[7,9],ext:7,extern:[9,11],extra:[11,17],f_miss:[2,11],fail:[2,5,8,9,16],fakedata:6,fals:[1,2,6,7,9],fam:[2,6,7,13,14,16,17],famili:11,far:11,fast:[9,17],fasta:7,father:11,featur:[3,9],femal:[2,11],feolo:2,few:[4,5,6,9],fid1:11,fid2:11,fid:11,field:[7,11],fifth:11,fig:10,fight:9,file1:9,file2:9,file:[1,2,3,4,5,6,10,12,13,14,15,16,17],file_nam:[1,2,5],file_pattern:[1,5],file_typ:7,filenam:[1,2,5],filenotfound:9,filenotfounderror:6,files_for_lab:13,filetyp:11,fill:[1,5,7],filter:[1,2,5,7,10],find:[4,9,11],findabl:4,first:[2,4,8,9,11,12,16,17],fix:9,fixtur:6,flag:[6,7,9,16,17],flake8:4,flight:[7,16],float32:2,flow:5,fna:7,folder:[3,5,6,8,9],follow:[2,3,4,5,6,9,10,11],forc:8,forg:[4,8],form:9,format:[2,4,5,7,11,14],formatt:4,forward:2,forward_strand_annot:2,found:[2,9,16],four:[11,14],fourth:11,frame:[2,5],framework:5,frequenc:[2,7,11,12],fresh:4,from:[1,2,3,4,5,6,7,9,11,12,16,17],full:[1,2,4,5,6,9,10,11],full_sample_sheet:6,fulltext:2,fundament:9,fwd:7,gatk:2,gc10:2,gc50:2,gca_000001405:7,gcc:4,gencal:[2,7],gender:[2,11],gener:[2,3,6,7,9,10,13,16,17],genet:11,geno:[2,7,11],genom:[2,7,9,12,16],genome_build:[7,9],genomebuild:[9,11],genomestudio:11,genotyp:[2,7,10,17],genotypecal:2,get:[2,3,5,9,10],get_autocall_d:2,get_autocall_vers:2,get_ballele_freq:2,get_base_cal:2,get_base_calls_forward_strand:2,get_base_calls_gener:2,get_base_calls_plus_strand:2,get_call_r:2,get_cluster_fil:2,get_control_x_intens:2,get_control_y_intens:2,get_gc10:2,get_gc50:2,get_gend:2,get_genotyp:2,get_genotype_scor:2,get_imaging_d:2,get_logr_dev:2,get_logr_ratio:2,get_normalization_transform:2,get_normalized_intens:2,get_num_cal:2,get_num_intensity_onli:2,get_num_no_cal:2,get_num_snp:2,get_percentiles_i:2,get_percentiles_x:2,get_ploidi:2,get_ploidy_typ:2,get_raw_x_intens:2,get_raw_y_intens:2,get_record_problem:2,get_sample_nam:2,get_sample_pl:2,get_sample_wel:2,get_scanner_data:2,get_slide_identifi:2,get_snp_manifest:2,gether:9,getting_start:3,giant:5,git:[3,4,6],github:[2,3,4,5,6,8,9,12,13,14,15,16,17],give:[1,5],given:[1,5,6,7,11],goal:5,going:[4,5],graf:[0,5,16,17],graf_ancestri:16,grch_version:6,gre:8,green:7,group:[1,5,9,16],group_bi:7,group_by_subject_id:[9,16],gsa_lab_qc:9,gsamd:[7,11],gtc2plink:7,gtc:[2,5,7,9,12,14,16],gtc_pattern:[7,9,12,14],gwa:[7,11],gwas_primaryqc:7,gwas_qc_log:9,gwasqcpipelin:[1,3,6,7,9,12,13,14,15,16,17],gwasqcpipeline_submiss:[5,9],gwasqcpipeline_v1:8,gwasqcpiplin:9,had:[5,9,16],half:2,hand:11,handl:[2,5],haploid:[2,11],haploview:11,hardi:2,harvard:2,has:[2,4,5,6,7,11,16],hat:7,have:[1,2,3,4,5,6,7,8,9,11,12,15,16,17],head:[1,5],header:[7,11],help:[4,8,9],helper:[2,3,5,6,10],here:[2,3,5,7,9,10,11,16,17],het:[2,7,11,17],hetconc:2,heterozyg:[2,11],heterozygos:[7,16],heterozygot:[2,11],hethet:[2,11],hg37:[7,9],hg38:[7,9],hg_match:2,hg_miss:2,hgmr:2,hide:9,high:[11,16],highli:16,hom:[2,11],home:[1,3,5],homhom:[2,11],homibs0:2,homo:11,homozyg:[2,11],homozygot:[2,17],hour:[7,9],how:[2,5,7,9,11],howev:[3,4,5,6,7],hsph:2,html:[2,3,4,5,11],http:[2,4,6,8,9,11,12,13,14,15,16,17],human:[7,9],hwe:[2,7,17],hwp:[7,13],ibc:2,ibd1:2,ibd1seg:2,ibd2:[2,16],ibd2seg:2,ibd:[2,7,11,16,17],ibd_pi_hat_max:[7,16,17],ibd_pi_hat_min:[7,16,17],ibs0:[2,11],ibs1:[2,11],ibs2:[2,11],id1:2,id2:2,idat:[5,7,9,15,16],idat_pattern:[7,9,14,15],idatintens:16,ident:2,identifi:[2,7,11,16,17],identifier_sex:11,identifil:[13,16],identifiler_need:16,identifiler_reason:16,ignor:[4,7,9,11],ignore_vers:2,iid1:[2,11],iid2:[2,11],iid:[2,11],illumina:[0,5,7,9],illumina_cluster_fil:7,illumina_manifest_fil:[7,12],illuminaadpcfilewrit:2,illuminabeadarrayfil:[2,5],illuminaio:5,ilmn_id:2,ilmnid:[2,11],ilmnstrand:11,imag:[2,3],imiss:2,implement:[3,9],imput:[2,11],inbreed:[2,11],includ:[2,5,6,7,9,11,16],inconveni:4,increas:9,indel:2,indep:7,independ:10,index:[2,3,7],indic:[7,9,11,17],individu:[2,7,11,12,15],infer:11,infinium:[2,11],info:[2,4],inform:[2,6,7,8,11],ini:7,init:4,initi:[8,16],input:[2,5,9,11],instal:[3,5,9,10],instanc:1,instead:[9,11],instruct:[4,8],int16:2,int32:2,integ:[2,7,11],integr:6,intens:[2,7,10,11,16],intensity_check:15,intensity_onli:11,intensity_threshold:[7,16],interact:[5,8,9],interfac:2,intern:[1,5,6,7,9],internalqcknown:16,interpret:[4,11],is_call_rate_filt:16,is_contamin:16,is_cr1_filt:16,is_cr2_filt:16,is_discordant_repl:16,is_ge_concord:16,is_ge_pi_hat:16,is_internal_control:[9,16],is_missing_gtc:[9,16],is_missing_idat:[9,16],is_sample_exclus:[9,16],is_sex_discord:16,is_subject_repres:[16,17],is_user_exclus:[7,9,16],is_write_complet:2,isn:11,isol:6,isort:4,issu:[4,9],item:2,iter:2,its:[6,9,16],itself:4,javadoc:2,jin:2,jinja2:5,job:[5,7,9],job_id:9,jobscript:5,json:3,jun:11,just:[2,5,7,9,11,13],keep:[4,5,7,9],kei:[5,9],keyr:4,keyword:3,kill:9,king:[0,5,16],knownrepl:[13,16],kwarg:[2,9],lab:16,label:11,laboratori:7,languag:3,languageserv:3,larg:[2,6,9,11],last:[2,11],latest:8,ld_prune_r2:[7,16,17],least:6,left:11,legaci:5,length:2,let:[4,10],level:[2,5,9,10,11,16],librari:[2,10],light:9,like:[3,7],likelihood:11,lim:[5,9],limit:[9,11],lims_individual_id:11,lims_output_dir:[7,13],lims_specimen_id:11,lims_upload:[7,13],limssample_id:11,limsupload:[7,13],line:[4,9,11],linear:11,link:[3,7,11,14],linter:[3,4],linux:[4,8],list:[1,2,5,7,8,9,11,16,17],lit:2,live:5,lln0:11,lln:11,lmiss:2,load:[8,11],load_config:[1,5],local:[4,6],local_mem_mb:9,local_task:9,locat:[2,3,4,5,7,8,9,11,13],loci:2,locu:2,locusentri:2,log:[9,11],logic:5,logist:7,logr:2,longer:9,lookup:2,loop:5,lot:5,low:[11,16],lscratch:8,made:[6,11,17],maf:[2,7,11,17],maf_for_hw:[7,17],maf_for_ibd:[1,5,7,16,17],magic:11,mai:[4,5,7,9],main:[3,4,5,7,9],mainten:5,major:[2,5,11,12,13,14,15,16,17],make:[2,3,4,5,6,9,16],make_cgr_sample_sheet:6,make_config:6,makefil:3,male:[2,11],mamba:[4,8],manag:[0,2,4,10],mani:[5,6],manichaikul:2,manifest:[2,5,9],manifest_nam:2,manual:2,map:[2,7,14],map_info:2,mapinfo:11,mark:[6,9],markdown:[3,5],marker:[2,7,11],marker_:11,markup:3,match:[2,3,4,7,11,16],materi:3,max:7,max_time_hr:7,maximum:7,mayb:9,mean:[5,11],median:[7,15,16],median_idat_intens:15,mem:[8,9],memori:9,merg:14,merlin:11,mess:2,metadata:[5,7,11,17],method:[1,2,5,16,17],metric:[10,16],microarrai:7,microsoft:2,mimic:6,min:7,mind:7,miniconda3:[3,4,8],miniconda:[4,8],minim:[5,11],minimum:[7,17],minimum_pop_subject:[7,17],minor:[2,7,11],minu:11,minut:4,mismatch:2,miss:[2,9,16,17],miss_pheno:[2,11],missing:16,mistak:4,mix:[7,11,16],mixtur:[11,16],mkdir:8,mode:[2,9],model:[5,6,7],modern:4,modif:16,modifi:9,modul:[1,5,6,8],module_dir:[1,5],moment:[2,17],mondai:2,monomorph:11,month:11,more:[3,4,5,6,7,9,11,12,13,14,15,16,17],morgan:[2,11],most:[5,7,11],mostli:[5,9,13,16],mother:11,move:5,msdn:2,mtdna:11,much:6,multipl:[7,10,16],must:[1,2,3,5,6,7,9,11],mutual:7,my_project_nam:9,mychaleckyj:2,myfold:8,mypi:4,n_clst:11,n_geno:[2,11],n_miss:[2,11],n_nm:[2,17],n_snp:2,name:[1,2,5,6,7,9,11,16,17],namesapac:7,namespac:[5,7],nan:16,ncbi:2,nchrob:11,nci:[4,6,8,12,13,14,15,16,17],ndarrai:2,necessari:[4,5,9],need:[3,4,5,6,7,8,9,11,16],neg:11,never:[3,5],newlin:2,next:[4,8,9,11],nice:5,nicer:5,nih:6,ninth:2,nnumber:4,non:[2,9,11,17],non_existing_fil:6,none:[1,2,4,5,6,8,9],nonetyp:2,nonexist:11,nonmiss:[2,11],nonzero:2,normal:[2,9,11],normalization_lookup:2,normalizationtransform:2,not_major_chrom:2,note:[2,9],notemp:9,notic:6,now:4,nsertion:2,nucleotid:2,num:7,num_analytic_exclus:16,num_loci:2,num_sampl:[7,9],num_samples_per_subject:[9,16],num_snp:[7,9],number:[2,4,5,6,7,9,10,11,16,17],numer:11,numpi:2,o_het:2,o_hom:[2,17],object:[2,5,7],obligatori:[2,11],observ:[2,11,17],obtain:2,offici:[4,11],offspr:2,often:[5,6],older:11,onc:[4,5,6,8],one:[2,4,5,6,7,9,11],onli:[1,2,4,5,6,7,9,11,17],open:2,optim:9,option:[1,5,6,7,9,10,11,12,13,14,15,16,17],orchestr:[6,9],order:[2,4,11],org:[2,4,11],origin:11,other:[2,3,5,6,7,9,11,16],otherwis:[2,3,8],our:[4,5,6,7,9,16],out:[2,4,5,7,17],output:[2,5,7,9,10,11,12,13,14,15,16,17],output_pattern:[7,13],overrid:9,overridden:2,own:[5,6,9,16],p_valu:2,packag:[2,3,4,5,10],page:[3,5,7,11],pair:[2,11,17],pairwis:[2,7,16],panda:[1,2,5],pandoc:5,panel:[11,17],param:[2,7],paramet:[1,2,5,6,11],parent1:11,parent2:11,parent:2,pars:[2,5,16],parser:[0,3,6,17],part:[5,6,7],parti:[5,7],particularli:9,partit:[7,9],partition_nam:9,pass:[2,5,6,7,9,16],path:[1,2,3,4,5,6,7,9,14],pathlib:[1,2,5,6],pathlik:[2,5],pattern:[1,5,7,9],payload:5,pc10:[2,17],pc1:[2,11,17],pc2:[2,11,17],pc3:17,pc4:17,pc5:17,pc6:17,pc7:17,pc8:17,pc9:17,pca:7,pca_plot:17,ped:[2,7,14],pedigre:11,pedsex:[2,11],peopl:2,per:[1,9,11,14,16,17],percent:[11,16],percentil:2,perfer:7,perform:7,permit:11,phase3_shapeit2_mvncall_integrated_v5:7,phase:9,phe:[2,11],phenotyp:[2,7,11,16,17],physic:11,pi_hat:[2,11,16,17],pi_hat_threshold:[7,16,17],pi_study_id:11,pi_subject_id:11,picard:[2,11],pictur:10,piec:6,pip:8,pipelin:[7,10,11],pipeline_vers:7,place:7,placehold:9,plain:5,plate:2,pleas:[2,9],plink2:[5,11],plink:[0,5,7,16,17],plink_merg:5,plo:2,ploidi:[2,11],plot:17,plot_autosomal_heterozygos:7,plu:[2,7,11],plugin:[3,5],png:[16,17],point:[4,5,6,7,10],pool:2,pop:[7,17],popgroup:11,popul:[7,9,10,11,16],population_level:17,population_level_ibd:7,population_qc:17,pos:2,posit:11,posixpath:[1,5],possibl:[3,6,9,16],potenti:[2,11],poulat:17,power:6,ppc:[2,11],pre:[7,16],pre_flight:5,predict:[2,11,16],predicted_sex:16,prefer:3,prefix:[3,7,12],prepend:[1,5],present:[7,9,16],pretti:5,preview:3,previou:10,price:2,princip:[2,17],principl:5,probabl:[2,5,9],probe:[2,11],problem:[2,4,5,9,11,16],problemat:11,process:9,prod:7,produc:[2,17],product:[6,8,9],profil:9,program:[5,11],project:[4,5,7,9,11],project_nam:[7,9],prompt:9,proper:9,properli:11,properti:[1,5,7],propibd:2,proport:[2,11,16],provid:[2,5,6,7,8,9,10,14,16],prune:[7,17],psutil:4,publicli:6,pull:[5,7,12,17],purpos:9,push:[3,4],pwd:4,py38_4:4,py39_4:8,py3:[4,8],pydant:[3,5],pyproject:4,pysam:4,pytest:[1,4,6],python:[1,3,4,5,8,10],python_keyring_backend:4,qc_exclus:5,qc_family_id:17,qc_report:[5,13],qc_report_data_dictionari:1,qc_v:9,qsub:9,quantit:11,queri:[1,5],queue:9,quickli:[2,6],qwasqcpipelin:4,rais:[6,9],ran:16,rang:2,rate:[2,7,11],ratio:[2,11],raw:2,reach:9,read:[2,9],read_byt:2,read_char:2,read_float:2,read_genom:2,read_het:[2,17],read_hw:2,read_imiss:2,read_int:2,read_lmiss:2,read_rel:2,read_related:2,read_scanner_data:2,read_sexcheck:2,read_str:2,read_ushort:2,readabl:9,readdata:6,reader:11,readthedoc:9,real_data:6,realdata:6,realiti:7,realli:[4,11],reason:4,recent:9,recomend:4,recommend:8,record:2,recurs:[4,6],red:7,ref_strand:2,ref_strand_annot:2,refer:[2,3,5,9,10,17],referenc:[5,7],reference_fasta:7,reference_fil:[5,7,12,16],refstrand:[2,11],regist:5,regress:7,reinstal:8,rel:[2,7,17],relat:[2,7,16,17],related:[10,16],related_subject:17,related_subjects_to_remov:17,relationship:11,releas:[3,4,8],reli:16,remain:11,remov:[2,16],remove_contam:[7,16],remove_rep_discord:[7,16],renam:17,rep:7,repeatedli:5,replac:6,replic:[6,7,11],replicate_id:[9,16],repo:[4,8],report:[1,2,6,7,11,13],report_strand:2,repositori:[4,6],repres:[2,7,11,16],requir:[4,5,6,7,8,9,10,12,15],resid:5,resolv:5,resourc:9,respect:[2,5,11],rest:11,restart:4,restructuredtext:3,result:[7,11,12,16,17],review:9,rich:2,right:11,road:6,robust:2,robustli:2,root:[1,3,5,7],round:11,row:[2,7,11],rsid:16,rst:3,rstcheck:[3,4],rsync:6,rule:[1,5,7,9],run:[3,4,5,6,7,8,10,11,16,17],runner:[1,5],runtim:[4,5],s0002:2,safe:11,sale:2,same:[7,11],samp0001:7,samp0002:7,samp0003:7,samp0004:7,sampl:[1,2,5,6,9,10,12,13,14,15,17],sample0001:7,sample_call_rate_1:[7,16],sample_call_rate_2:[7,16],sample_concordance_plink:7,sample_group:11,sample_id:[1,5,7,9,11,16,17],sample_ids_to_remov:7,sample_level:[7,12,14,15,16,17],sample_nam:11,sample_pl:11,sample_qc:[5,9,16,17],sample_qc_t:[16,17],sample_sheet:[5,7,9],sample_sheet_fil:[1,5,7],sample_typ:11,sample_wel:11,sampleusedforsubject:13,saniti:11,sankefil:5,sapien:11,sas_af:7,save:[5,6,7],sb001:[1,5],sb002:[1,5],sbatch:9,scan:2,scanner:2,scannerdata:2,schedular:9,score:[2,7,11,12],scratch:8,script:[1,5,7,9,16,17],scripts_dir:[1,5],search:[2,9,16],second:[2,9,11,16],section:[5,7],see:[1,2,4,5,6,7,8,9,12,13,14,15,16,17],segment:2,select:7,self:2,send:9,sentrixbarcode_a:11,sentrixposition_a:11,separ:[2,4,11,17],sequenc:11,seri:[5,16],server:6,session:[1,8],set:[1,2,3,5,6,7,9,10,11,13,16,17],setup:[3,4],seventh:11,sever:[5,10],sex:[2,7,9,16],sex_chr_includ:7,sex_dtyp:16,sexcheck:2,sge:9,share:[2,5,6,11],she:11,sheet:[1,5,6,9,16],sherri:2,should:[2,4,6,7,8],show:4,shrinkag:11,sib:2,sibl:2,side:9,similar:3,simpl:[6,11],simpli:[6,17],simplifi:9,sinc:11,singl:[2,3,7,17],sinteract:8,site:7,six:[10,11],sixth:11,size:[9,11,17],skip:[4,6,9],slow:4,slurm:[7,9],slurm_partit:[7,9],small:[4,6,11],smartpca:2,smk:[5,12,13,14,15,16,17],snakefil:[1,5,9],snakefmt:4,snakemak:[1,5,10],snp1:11,snp2:11,snp3:11,snp:[2,7],snp_arrai:7,snp_call_rate_1:[7,16],snp_call_rate_2:[7,16],snp_cr1:16,snp_cr1_remov:16,snp_cr2:16,snp_cr2_remov:16,snp_id:11,snp_qc:16,snp_qc_tabl:16,snpsex:[2,11],snpweight:16,softwar:[2,5,6,11],software_param:[1,5,7,12,16,17],solv:5,some:[5,6,9],someth:[4,5],sometim:[4,9],sort:2,sourc:[1,2,3,4,6,8,9,10,11],source_strand:2,sourceseq:11,sourcestrand:[2,11],sourcevers:11,sp001:[1,5],sp002:[1,5],sp003:[1,5],sp00:[1,5],space:5,spec:11,speci:11,special:[6,7],specif:[6,7,9,11],specifi:[5,7,9,11],spell:4,sphinx:3,sphinxbuildpath:3,split:[10,17],sqrt:17,squar:7,sr0001:7,sr_subject_id:11,src:[1,5,6,12,13,14,15,16,17],src_dir:[1,5],ssl:4,stabl:9,stage:4,standard:[9,11],start:[3,7,8,10],state:[2,6],statist:11,statu:[2,5,7,9,11,16,17],step:[4,6,9,10],still:11,store:[3,4,5,6,11],str:[1,2,5,9],strand:[2,7,11],strand_annot:2,strategi:9,string:[2,7,16,17],structur:[6,7,9,10],studi:2,studysampleknown:16,style:3,sub0001:7,sub0002:7,sub0003:7,sub:[1,5,7,9,10],sub_workflow:[1,3,5,12,13,14,15,16,17],subclass:2,subdirectori:6,subfold:5,subject:[2,5,7,9,10,13,16],subject_dropped_from_studi:16,subject_id:[1,5,7,17],subject_id_column:[7,9],subject_level:17,subject_qc:[5,9,17],subjects_maf:17,submiss:[5,9],submission_cmd:9,submit:[5,7],submodul:[4,6],subset:16,substitut:5,subworkflow:[1,5,9],subworkflow_dir:[1,5],success:5,successfulli:2,suggest:[5,7,8,9],suit:[4,10],summar:16,summari:10,summary_stat:[5,16],support:[2,7,11],supported_vers:2,sure:[4,5,9,11,16],svalid:9,symbol:[7,14],symlink:6,sync:6,synthet:6,system:[4,5,6,7,9,10],tabl:[5,7,11,12,13,15,16],tag:4,take:[4,9,10,14,16],tar:4,target:[6,9],task:[5,9],tbi:[7,9],techniqu:6,tell:[4,5,6,9],templat:[1,5],template_dir:[1,5],temporari:9,tend:4,tenth:2,termin:4,terminolog:5,test:[2,4,5,7,9,10,11],test_:6,test_addit:6,test_config:6,test_data:6,test_data_path:6,test_data_serv:6,test_data_us:6,test_snakemake_job_group:6,text:[5,11],than:11,thead:9,thei:[4,6,7,9,16],them:[4,5,9,11,14,17],thi:[1,2,3,4,5,6,7,9,10,11,12,13,15,16,17],thing:[1,4,5],third:[2,11],those:[2,6],though:7,thought:4,thousand:7,thousand_genom:5,thousand_genome_:7,thousand_genome_tbi:[7,12,16],thousand_genome_vcf:[7,12,16],thousand_genomes_snp_id:16,thread:9,three:[2,4,5,7,9,11,16],threshold:[7,16,17],through:7,throughout:[3,5,11],thu:11,tidi:9,time:[4,7,9],time_hr:9,time_start:[7,13],titl:7,tmp:8,togeth:[14,17],toml:4,tool:[4,5,7,9,10],toolset:2,top:[2,5,10,11],topgenomicseq:11,total:[2,16],track:4,tradition:5,trait:11,transform:[2,11],translat:9,tri:[5,9],trigger:[3,7],tupl:2,twin:2,two:[5,7,9,11,16,17],txt:[1,5,6,12,16],type:[2,3,4,5,7,10],typer:5,typic:9,uint32:[2,16],unaff:[2,11],understand:11,union:2,uniqu:7,unit:6,unknown:[2,7,11,16,17],unknown_annot:2,unknownrepl:[13,16],unless:[6,7,9],unlik:2,unrel:2,unsur:2,unus:9,updat:[3,4,7,8,9],upload:7,upon:2,url:3,use:[1,2,3,4,5,6,7,8,9,11,16,17],used:[2,3,5,6,7,9,11,16,17],useful:9,user:[1,3,5,6,9,10],user_config:[1,5],user_fil:[5,7,9,12,13,14,15],uses:[3,4,5,6],ushort:2,using:[2,3,4,5,6,7,9,12,14,16],usual:[2,11],util:[5,9,10],valid:[2,6,9,11],valu:[5,6,7,9,11],variabl:6,variant:[2,11],variou:[5,7,10,11,16],vcf:[5,7,9,12,16],vector:2,venv:4,veri:[5,6,9],verifi:[2,4],verifyidintens:[5,12],version:[2,4,5,6,7,8,9,11],virginia:2,virtualenv:4,vvv:4,wai:[3,5,6,7],walk:7,walltim:9,want:[4,5,7,9],wc9c:2,web:5,websit:[2,4,7,8],weight:[4,9,11],weinberg:[2,11],well:[2,9,11],were:[2,16],wget:[4,8],wgs:7,what:[5,7,9,11],when:[2,3,5,7,9,11],where:[6,7,9,11,16,17],whether:2,which:[2,3,5,7,9,10,11,16],whl:[4,8],why:5,wide:2,wildcard:[5,7],within:[2,11],without:[5,7],work:[1,4,7,9,10],workflow:[1,3,4,6,8,10],workflow_dir:[1,5],workflow_param:[5,7,9,13,16,17],working_dir:6,workspacefold:3,would:[5,6,7,9,16],wrapper:9,write:[2,5],written:[2,3],wrote:5,www:[2,11],x86_64:[4,8],x_inbreeding_coeffici:16,xlsx:[1,13],yaml:[4,5],ycount:11,year:11,yield:2,yml:[1,5,9,12,13,14,15,16,17],you:[2,3,4,5,6,7,8,9,11,12,15,16],your:[3,4,5,6,7,8,9,11],yri:11,yzxa6408:2,zcat:7,zero:[2,11],zip:[1,5,13]},titles:["API Reference","Configuration API","Parsers","Working with Documentation","Development Environment Set-up","Working with Source Code","Working with Test Suite","Configuration Files","Installing GwasQcPipeline","Running the Pipeline","GwasQcPipeline Documentation","File Type Reference","Contamination Sub-workflow","Delivery Sub-workflow","Entry Points Sub-workflow","Intensity check Sub-workflow","Sample QC Sub-workflow","Subject QC Sub-workflow"],titleterms:{IDs:7,The:7,abf:11,adpc:11,analysi:17,ancestri:16,api:[0,1],autosom:17,bed:11,bim:11,bin:11,biowulf:8,bpm:11,build:3,call:[11,16],ccad2:8,ccad:8,cgem:8,cgr:[7,9],check:[9,11,15],cli:5,cluster:[5,9],code:5,command:5,commit:4,concord:16,conda:4,config:[7,9],configur:[1,5,7,9],consist:4,contam:11,contamin:[12,16],creat:[4,9],csv:11,data:6,deliveri:13,depend:4,develop:4,document:[3,10],download:4,eigensoft:[2,11],eigenstratgeno:11,entri:14,environ:4,exampl:7,fake:6,fam:11,file:[7,9,11],filter:16,fingerprint:2,flight:9,freq:11,frq:11,full:7,gener:5,genet:2,genom:11,genotyp:11,graf:2,gtc:11,gwasqcpipelin:[4,5,8,10],hardi:[11,17],heterozygos:17,hook:4,hwe:11,idat:11,illumina:[2,11],imiss:11,ind:11,indep:11,infer:2,instal:[4,8],intens:15,interfac:5,intern:[16,17],king:2,kinship:2,level:[7,17],lim:7,line:5,lmiss:11,local:[3,9],make:11,manag:[1,5],manifest:[7,11],map:11,miss:11,nih:[5,8],other:8,out:11,pairwis:11,paramet:7,parser:[2,5],pca:17,ped:11,pipelin:9,plink:[2,11],poetri:4,point:14,popul:17,pre:[4,9],profil:5,prune:11,publish:3,rate:16,real:6,refer:[0,7,11],related:[2,17],relationship:2,remov:7,replic:16,report:[5,16,17],repres:17,run:9,sampl:[7,11,16],select:17,set:4,sex:11,sexcheck:11,sheet:[7,11],snakemak:9,snp:[11,16],snpweight:11,snpwt:11,softwar:7,sourc:5,structur:17,sub:[12,13,14,15,16,17],subject:17,submit:9,suit:6,summari:[16,17],system:8,tabl:[2,17],test:6,top:7,txt:11,type:11,user:7,valid:5,valu:2,verifyidintens:11,virtual:4,weinberg:17,work:[3,5,6],workflow:[5,7,9,12,13,14,15,16,17],yml:7}})