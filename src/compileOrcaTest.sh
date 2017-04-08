filepath=/home/max/Projects/OrcaVMD/water.out
if [ "$1" ]; then

if [ $1 == "f" ]; then
	echo "reading fake input"
	filepath=/home/max/Projects/OrcaVMD/fake.out
fi
fi
gcc main.c orcaplugin.c -I../../include -o plugintest -lm
./plugintest orca $filepath
