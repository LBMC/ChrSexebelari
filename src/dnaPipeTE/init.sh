mkdir bin
cd bin

# install gnu parallel
sudo apt-get install parallel
cp $(which parallel) parallel

# install trf
curl -k -L http://tandem.bu.edu/trf/downloads/trf409.linux64 -o trf
 trf

# RM blast
curl -k -L  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28/ncbi-rmblastn-2.2.28-x64-linux.tar.gz -o ncbi-rmblastn-2.2.28-x64-linux.tar.gz
tar -xvf ncbi-rmblastn-2.2.28-x64-linux.tar.gz

# blast
curl -k -L ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/ncbi-blast-2.2.28+-x64-linux.tar.gz -o ncbi-blast-2.2.28+-x64-linux.tar.gz
tar -xvf ncbi-blast-2.2.28+-x64-linux.tar.gz

# RM
curl -k -L http://repeatmasker.org/RepeatMasker-open-4-0-6.tar.gz -o RepeatMasker-open-4-0-6.tar.gz
tar -xvf RepeatMasker-open-4-0-6.tar.gz
curl -k -L http://www.ens-lyon.fr/LBMC/intranet/fichiers/bioinfo/RepeatMaskerConfig.pm/at_download/file -o RepeatMasker/RepeatMaskerConfig.pm

# RM database
curl -k -L http://www.ens-lyon.fr/LBMC/intranet/fichiers/bioinfo/RM_Library.tar.gz/at_download/file -o RM_Library.tar.gz.tar.gz
tar -xvf RM_Library.tar.gz.tar.gz
rm -Rf RepeatMasker/Libraries
mv Libraries RepeatMasker/
var="#!"$(which perl)
echo $var | cat - RepeatMasker/util/queryTaxonomyDatabase.pl > temp && mv temp RepeatMasker/util/queryTaxonomyDatabase.pl
echo $var | cat - RepeatMasker/util/rmOutToGFF3.pl > temp && mv temp RepeatMasker/util/rmOutToGFF3.pl
echo $var | cat - RepeatMasker/util/rmToUCSCTables.pl > temp && mv temp RepeatMasker/util/rmToUCSCTables.pl
echo $var | cat - RepeatMasker/util/queryRepeatDatabase.pl > temp && mv temp RepeatMasker/util/queryRepeatDatabase.pl
echo $var | cat - RepeatMasker/RepeatMasker > temp && mv temp RepeatMasker/RepeatMasker
echo $var | cat - RepeatMasker/ProcessRepeats > temp && mv temp RepeatMasker/ProcessRepeats
echo $var | cat - RepeatMasker/RepeatProteinMask > temp && mv temp RepeatMasker/RepeatProteinMask
echo $var | cat - RepeatMasker/DateRepeats > temp && mv temp RepeatMasker/DateRepeats
echo $var | cat - RepeatMasker/DupMasker > temp && mv temp RepeatMasker/DupMasker

chmod +x RepeatMasker/util/queryTaxonomyDatabase.pl RepeatMasker/util/rmOutToGFF3.pl RepeatMasker/util/rmToUCSCTables.pl RepeatMasker/util/queryRepeatDatabase.pl RepeatMasker/RepeatMasker RepeatMasker/ProcessRepeats RepeatMasker/RepeatProteinMask RepeatMasker/DateRepeats RepeatMasker/DupMasker

# trinity
curl -k -L http://www.ens-lyon.fr/LBMC/intranet/fichiers/bioinfo/trinityrnaseq_r2013_08_14.tar.gz/at_download/file -o trinityrnaseq_r2013_08_14.tar.gz
tar -xvf trinityrnaseq_r2013_08_14.tar.gz
