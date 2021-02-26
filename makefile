git clone https://github.com/neurobin/shc.git
cd shc/
./autogen.sh
./configure
make
cd ..
./shc/src/shc -f ./bash_scripts/genome_download.sh -o genome-download
./shc/src/shc -f ./bash_scripts/plot-variants.sh -o plot-variants
./shc/src/shc -f ./bash_scripts/variants2genes.sh -o variants2genes
mkdir bin
mv genome-download plot-variants variants2genes ./bin/
echo "bin folder containing executable binaries are made"
echo ""
echo "make done"
