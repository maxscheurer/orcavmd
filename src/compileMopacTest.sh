filepath=$(pwd)/mopac_testfiles/test.out
if [ "$1" ]; then

if [ $1 == "f" ]; then
	echo "reading fake input"
	filepath=$(pwd)/mopac_testfiles/fake.out
elif [ $1 == "s" ]; then
	echo "reading wrong version file"
	filepath=$(pwd)/mopac_testfiles/unsupp.out
else
	filepath=$1
fi
fi
echo $filepath

set -e
g++ main.C -g -I../../include -c -o main.o
# g++ Matrix.cpp -g -std=c++11 -c -o Matrix.o
g++ mopacplugin.C -g -std=c++11 -isystem../../include -c -o plugin.o
gcc plugin.o main.o -lstdc++ -lm -o plugintest
./plugintest mopac $filepath

#make LINUXAMD64 -j9; sudo cp compile/lib_LINUXAMD64/molfile/orcaplugin.so /usr/local/lib/vmd/plugins/LINUXAMD64/molfile/
