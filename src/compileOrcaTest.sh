filepath=$(pwd)/orca_testfiles/water.out
if [ "$1" ]; then

if [ $1 == "f" ]; then
	echo "reading fake input"
	filepath=$(pwd)/orca_testfiles/fake.out
elif [ $1 == "s" ]; then
	echo "reading wrong version file"
	filepath=$(pwd)/orca_testfiles/unsupp.out
else
	filepath=$1
fi
fi
echo $filepath
gcc main.c orcaplugin.c -I../../include -o plugintest -lm
./plugintest orca $filepath

#make LINUXAMD64 -j9; sudo cp compile/lib_LINUXAMD64/molfile/orcaplugin.so /usr/local/lib/vmd/plugins/LINUXAMD64/molfile/
