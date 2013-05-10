#!/bin/bash 

HOME_DIR=$PWD
TEST_DIRECTORY=${1}
TEST_FLAG=${2}

cd ../current/$TEST_DIRECTORY

COUNT=0
for x in *
do
	case $x in 
		*5_*.dot) 
			COUNT=$(( $COUNT + 1 ))
			NAMEFILE=$( printf 'stage5_%02d.pdf' $COUNT )
			echo $NAMEFILE
			if [ $TEST_FLAG -eq 1 ]; then
				dot -Tpdf $x -o $NAMEFILE
			fi

			if [ $TEST_FLAG -eq 2 ]; then 
				sfdp -Tpdf $x -o $NAMEFILE				
			fi

			if [ $TEST_FLAG -eq 3 ]; then
				sfdp -Goverlap=prism -Tpdf $x -o $NAMEFILE
			fi
		;;
	esac
done 


