for i in {90..92..2}
do
	mpirun -np 2 ./main $i 20
	mv './output.csv' "./output_$i.csv"
    cp ./output_$i.csv ./outputfile
	rm ./output_$i.csv	
done
