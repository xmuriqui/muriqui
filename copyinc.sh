#!/bin/bash


incDirs="codbb codminlpp codspm codopt iquad coddctools coddctools/pugixml codnumcomp"

curDir=`pwd`

#echo "current: " $curDir

mkdir includes


for d in $incDirs; do
  if [ -d $d ]
  then
	echo $d;
	cd $d;
	
	for incFile in *.h*; do
		echo $incFile
		ln -s -f ../$d/$incFile  $curDir/includes/$incFile;
	done
	cd $curDir;
  fi
done


cd $curDir

#copying the include files in the current directory
for incFile in *.h*; do
  echo $incFile
  ln -s -f ../$incFile  $curDir/includes/$incFile;
done


