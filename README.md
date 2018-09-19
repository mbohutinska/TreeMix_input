# TreeMix_input
Scripts to prepare TreeMix input from diploid or polyploid vcf input  and visualize the output


### Documentation converting vcf to TreeMix input, based on python script of Christian Sailer (Christian.Sailer@jic.ac.uk) and R script of Filip Kolar (filip.kolar@gmail.com)

### 1. to extract variable sites from foufold.filtered:
qsub filter4dVariable.sh 

### 2. run .splitVCFsTreeMix function ScanTools_ProtEvol  - MetaCentrum 
test.splitVCFsTreeMix(vcf_dir="arenosa_fourfold.variable", pops=['BAB','SUB','TRD'], mem="16", time_scratch='00:10:00', ncpu="5", scratch_path="$SCRATCHDIR",min_dp="8",mffg="0.2", overwrite=True, scratch_gb="8", keep_intermediates=False, use_scratch=True, print1=False)

### 3. copy .splitVCFsTreeMix to local folder, run:
python3 conversionTreemixMajda.py -i "Tatry/" -o "Tatry/"

### 4. run treemix
treemix -i Tatry/treemix_input.table.gz -o Tatry/outstem
build ML tree of the populations
 -root = set the outgroup, e.g. named CRO
 -k = set how many successive SNPs are in LD (here I assume just every one is independent, after previous filtering, thus setting -k 1)
 -m allow for m migration events
treemix -i DIPLOIDS.treemix.frq.gz -root CRO -k 1 -m 2 -o outstem

treemix -i RADtetDip/treemix_input.table.gz -root OUT -o RAD/RADtetDip_mig0_boot0 
treemix -i RADtetDip/treemix_input.table.gz -root OUT -m 1 -o RAD/RADtetDip_mig1_boot0 
treemix -i RADtetDip/treemix_input.table.gz -root OUT -m 2 -o RAD/RADtetDip_mig2_boot0 
treemix -i RADtetDip/treemix_input.table.gz -root OUT -m 3 -o RAD/RADtetDip_mig3_boot0 
treemix -i RADtetDip/treemix_input.table.gz -root OUT -m 4 -o RAD/RADtetDip_mig4_boot0 

###bootstrap over 100 replicates
##creates a file BOOTSTRAP.tre in the bootstrapping folder, could be opened e.g. in Figtree
#Without migration:
for i in {1..100}; do treemix -i treemix_input.table.gz -bootstrap -k 150 -o bootstrap/replicate$i; done
gunzip replic*treeout.gz
sumtrees.py --decimals=0 --percentages --output-tree-filepath=BOOTSTRAP.tre --target=replicate1.treeout replicate2.treeout replicate3.treeout replicate4.treeout replicate5.treeout replicate6.treeout replicate7.treeout replicate8.treeout replicate9.treeout replicate10.treeout replicate11.treeout replicate12.treeout replicate13.treeout replicate14.treeout replicate15.treeout replicate16.treeout replicate17.treeout replicate18.treeout replicate19.treeout replicate20.treeout replicate21.treeout replicate22.treeout replicate23.treeout replicate24.treeout replicate25.treeout replicate26.treeout replicate27.treeout replicate28.treeout replicate29.treeout replicate30.treeout replicate31.treeout replicate32.treeout replicate33.treeout replicate34.treeout replicate35.treeout replicate36.treeout replicate37.treeout replicate38.treeout replicate39.treeout replicate40.treeout replicate41.treeout replicate42.treeout replicate43.treeout replicate44.treeout replicate45.treeout replicate46.treeout replicate47.treeout replicate48.treeout replicate49.treeout replicate50.treeout replicate51.treeout replicate52.treeout replicate53.treeout replicate54.treeout replicate55.treeout replicate56.treeout replicate57.treeout replicate58.treeout replicate59.treeout replicate60.treeout replicate61.treeout replicate62.treeout replicate63.treeout replicate64.treeout replicate65.treeout replicate66.treeout replicate67.treeout replicate68.treeout replicate69.treeout replicate70.treeout replicate71.treeout replicate72.treeout replicate73.treeout replicate74.treeout replicate75.treeout replicate76.treeout replicate77.treeout replicate78.treeout replicate79.treeout replicate80.treeout replicate81.treeout replicate82.treeout replicate83.treeout replicate84.treeout replicate85.treeout replicate86.treeout replicate87.treeout replicate88.treeout replicate89.treeout replicate90.treeout replicate91.treeout replicate92.treeout replicate93.treeout replicate94.treeout replicate95.treeout replicate96.treeout replicate97.treeout replicate98.treeout replicate99.treeout replicate100.treeout
#With migration:
for i in {1..100}; do treemix -i treemix_input.table.gz -bootstrap -k 150 -o bootstrap/replicate$i; done
cd bootstrap
mkdir nomigr
for i in {1..100}; do head -n1 replicate$i.treeout > nomigr/replicate$i.treeout; done
cd nomigr
gunzip replic*treeout.gz
sumtrees.py --decimals=0 --percentages --output-tree-filepath=BOOTSTRAP.tre --target=replicate1.treeout replicate2.treeout replicate3.treeout replicate4.treeout replicate5.treeout replicate6.treeout replicate7.treeout replicate8.treeout replicate9.treeout replicate10.treeout replicate11.treeout replicate12.treeout replicate13.treeout replicate14.treeout replicate15.treeout replicate16.treeout replicate17.treeout replicate18.treeout replicate19.treeout replicate20.treeout replicate21.treeout replicate22.treeout replicate23.treeout replicate24.treeout replicate25.treeout replicate26.treeout replicate27.treeout replicate28.treeout replicate29.treeout replicate30.treeout replicate31.treeout replicate32.treeout replicate33.treeout replicate34.treeout replicate35.treeout replicate36.treeout replicate37.treeout replicate38.treeout replicate39.treeout replicate40.treeout replicate41.treeout replicate42.treeout replicate43.treeout replicate44.treeout replicate45.treeout replicate46.treeout replicate47.treeout replicate48.treeout replicate49.treeout replicate50.treeout replicate51.treeout replicate52.treeout replicate53.treeout replicate54.treeout replicate55.treeout replicate56.treeout replicate57.treeout replicate58.treeout replicate59.treeout replicate60.treeout replicate61.treeout replicate62.treeout replicate63.treeout replicate64.treeout replicate65.treeout replicate66.treeout replicate67.treeout replicate68.treeout replicate69.treeout replicate70.treeout replicate71.treeout replicate72.treeout replicate73.treeout replicate74.treeout replicate75.treeout replicate76.treeout replicate77.treeout replicate78.treeout replicate79.treeout replicate80.treeout replicate81.treeout replicate82.treeout replicate83.treeout replicate84.treeout replicate85.treeout replicate86.treeout replicate87.treeout replicate88.treeout replicate89.treeout replicate90.treeout replicate91.treeout replicate92.treeout replicate93.treeout replicate94.treeout replicate95.treeout replicate96.treeout replicate97.treeout replicate98.treeout replicate99.treeout replicate100.treeout

### 5. visualize the tree
treemix.R

### 6. run Figtree:
java -Xms64m -Xmx512m -jar /home/aa/Desktop/programs/FigTree_v1.4.3/lib/figtree.jar $*

