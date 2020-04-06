git clone https://github.com/neurobin/shc.git
cd shc/
./autogen.sh
./configure
make
cd ..
./shc/src/shc -f ./bash_scripts/sort_bam.sh -o sortBam
./shc/src/shc -f ./bash_scripts/genome_download.sh -o genomeDownload
./shc/src/shc -f ./bash_scripts/plot-variants.sh -o plotVariants
./shc/src/shc -f ./bash_scripts/variants2genes.sh -o variants2genes
mkdir bin
mv sortBam genomeDownload plotVariants variants2genes ./bin/
echo "bin folder containing executable binaries are made"
echo ""
echo "make done"