cp rmsd.txt rmsd.dat
cp matrix.txt matrix.dat
cp pop.txt pop.dat
cp matrix_element.txt matrix_element.dat

awk '{print $1","$3}' rmsd.dat > plot_pop_rmsd.csv
awk '{print $1 ","$5}' rmsd.dat > plot_co_rmsd.csv
sed -n '$p' pop.dat > temp

awk '{print $3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14}' temp > sample_input.dat

tail -12 matrix.dat > in_ma.dat



gfortran mc17_samples.f90 -o samples
./samples




gfortran p_corr.f90 -o p_corr
./p_corr

 
paste -d',' av_pop.csv  final_pop.csv >  plot_pop.csv
paste -d',' int_ma_initial.csv  int_ma_final.csv >  plot_int_ma.csv
paste -d',' int_ma_initial.csv  p_corr.txt >  plot_co_ma.csv
awk '{for(i=1;i<=NF;i++){printf "%s,", $i}; printf "\n"}' matrix_element.dat > plot_ma_element.csv
awk '{for(i=3;i<=NF;i++){printf "%s,", $i}; printf "\n"}' pop.dat > plot_pop_vs_time.csv



python python_subplot.py
