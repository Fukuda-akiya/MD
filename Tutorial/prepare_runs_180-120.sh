#!/bin/bash
#export CUDA_VISIBLE_DEVICES="0"

# Make Each Directory
for ((i=3; i<60; i+=3))
do
	mkdir -p min_$((180-i))
	mkdir -p min_$((180-i))/output
	mkdir -p equil_$((180-i))
	mkdir -p equil_$((180-i))/output
	mkdir -p samp_$((180-i))
	mkdir -p samp_$((180-i))/output
done

# Make Sampling File
for ((i=3; i<60; i+=3))
do
	echo "100 ps NPT production for $((180-i)) deg" > samp_$((180-i))/samp_$((180-i)).in
	echo " &cntrl" >> samp_$((180-i))/samp_$((180-i)).in  
	echo "  imin = 0," >> samp_$((180-i))/samp_$((180-i)).in  
	echo "  ntx = 5," >> samp_$((180-i))/samp_$((180-i)).in  
	echo "  irest = 1," >> samp_$((180-i))/samp_$((180-i)).in  
	echo "  ntpr = 10000," >> samp_$((180-i))/samp_$((180-i)).in  
	echo "  ntwr = 0," >> samp_$((180-i))/samp_$((180-i)).in  
	echo "  ntwx = 10000," >> samp_$((180-i))/samp_$((180-i)).in  
	echo "  ntf = 2," >> samp_$((180-i))/samp_$((180-i)).in  
	echo "  ntc = 2," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  cut = 8.0," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  ntb = 2," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  nstlim = 100000," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  dt = 0.001," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  temp0 = 300.0," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  ntt = 3," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  gamma_ln = 1," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  ntp = 1," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  pres0 = 1.0," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  taup = 5.0," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  nmropt = 1," >> samp_$((180-i))/samp_$((180-i)).in
	echo "  ioutfm = 1," >> samp_$((180-i))/samp_$((180-i)).in
	echo " &end" >> samp_$((180-i))/samp_$((180-i)).in
	echo " &wt" >> samp_$((180-i))/samp_$((180-i)).in
	echo "  type = "END"," >> samp_$((180-i))/samp_$((180-i)).in
	echo " &end" >> samp_$((180-i))/samp_$((180-i)).in
	echo "DISANG = samp_$((180-i))/disang.$((180-i))" >> samp_$((180-i))/samp_$((180-i)).in
	echo "DUMPAVE = dihedral_$((180-i)).dat" >> samp_$((180-i))/samp_$((180-i)).in
done

# Make Restrain File
for ((i=3; i<60; i+=3))
do
	echo "Harmonic restraints for $((180-i)) deg" > samp_$((180-i))/disang.$((180-i))
	echo " &rst" >> samp_$((180-i))/disang.$((180-i)) 
	echo "  iat = 9, 15, 17, 19" >> samp_$((180-i))/disang.$((180-i))
	echo "  r1 = $((0-i)), r2 = $((180-i)), r3 = $((180-i)), r4 = $((360-i))," >> samp_$((180-i))/disang.$((180-i))
	echo "  rk2 = 200.0, rk3 = 200.0," >> samp_$((180-i))/disang.$((180-i))
	echo "/" >> samp_$((180-i))/disang.$((180-i))
done

# Make Min Input File
for ((i=3; i<60; i+=3))
do
	echo "2000 step minimization for $((180-i)) deg" > min_$((180-i))/min_$((180-i)).in
	echo " &cntrl"  >> min_$((180-i))/min_$((180-i)).in
	echo "  imin = 1," >> min_$((180-i))/min_$((180-i)).in
	echo "  maxcyc = 2000," >> min_$((180-i))/min_$((180-i)).in
	echo "  ncyc = 500," >> min_$((180-i))/min_$((180-i)).in
	echo "  ntpr = 100," >> min_$((180-i))/min_$((180-i)).in
	echo "  ntwr = 1000," >> min_$((180-i))/min_$((180-i)).in
	echo "  ntf = 1," >> min_$((180-i))/min_$((180-i)).in
	echo "  ntc = 1," >> min_$((180-i))/min_$((180-i)).in
	echo "  cut = 8.0," >> min_$((180-i))/min_$((180-i)).in
	echo "  ntb = 1," >> min_$((180-i))/min_$((180-i)).in
	echo "  ntp = 0," >> min_$((180-i))/min_$((180-i)).in
	echo "  nmropt = 1," >> min_$((180-i))/min_$((180-i)).in
	echo " &end" >> min_$((180-i))/min_$((180-i)).in
	echo " &wt" >> min_$((180-i))/min_$((180-i)).in
	echo "  type = 'END'" >> min_$((180-i))/min_$((180-i)).in
	echo " &end" >> min_$((180-i))/min_$((180-i)).in
	echo "DISTANG = samp_$((180-i))/disang.$((180-i))" >> min_$((180-i))/min_$((180-i)).in
done

# Make Equil File
for ((i=3; i<60; i+=3))
do
	echo "50 ps NPT equilibration for $((180-i)) deg" > equil_$((180-i))/equil_$((180-i)).in
	echo " &cntrl" >> equil_$((180-i))/equil_$((180-i)).in
	echo "  imin = 0," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntx = 1," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  irest = 0," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntpr = 5000," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntwr = 50000," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntwx = 0," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntf = 2," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntc = 2," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  cut = 8.0," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntb = 2," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  nstlim = 50000," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  dt = 0.001," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  tempi = 0.0" >> equil_$((180-i))/equil_$((180-i)).in
	echo "  temp0 = 300.0," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntt = 3," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  gamma_ln = 1.0," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ntp = 1," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  pres0 = 1.0," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  taup = 5.0," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  nmropt = 1," >> equil_$((180-i))/equil_$((180-i)).in
	echo "  ioutfm = 1," >> equil_$((180-i))/equil_$((180-i)).in
	echo " &end" >> equil_$((180-i))/equil_$((180-i)).in
	echo " &wt" >> equil_$((180-i))/equil_$((180-i)).in
	echo "  type = 'END'" >> equil_$((180-i))/equil_$((180-i)).in
	echo " &end" >> equil_$((180-i))/equil_$((180-i)).in
	echo "DISTANG = samp_$((180-i))/disang.$((180-i))"  >> equil_$((180-i))/equil_$((180-i)).in
done

#################################	Conduct simulation	############################################
############################################################################################################
############################################################################################################
############################################################################################################
PMEMD=/home/fukuda/Software/amber20/bin/pmemd
for ((i=3; i<60; i+=3))
do      
	if [ $i -eq 3 ]; then
		$PMEMD -O -i min_$((180-i))/min_$((180-i)).in -o min_$((180-i))/output/min_$((180-i)).out \
			  -p tleap/ala_tri.prmtop \
		  	  -c equil/output/ala_tri_equil.rst7 -r min_$((180-i))/output/min$((180-i)).rst7
	   	
		$PMEMD -O -i equil_$((180-i))/equil_$((180-i)).in -o equil_$((180-i))/output/equil_$((180-i)).out \
			  -p tleap/ala_tri.prmtop \
			  -c min_$((180-i))/output/min$((180-i)).rst7 -r equil_$((180-i))/output/equil_$((180-i)).rst7

		$PMEMD -O -i samp_$((180-i))/samp_$((180-i)).in -o samp_$((180-i))/output/samp_$((180-i)).out \
			  -p tleap/ala_tri.prmtop \
			  -c equil_$((180-i))/output/equil_$((180-i)).rst7 -r samp_$((180-i))/output/samp_$((180-i)).rst7 \
			  -x samp_$((180-i))/output/samp_$((180-i)).nc

	else
		$PMEMD -O -i min_$((180-i))/min_$((180-i)).in -o min_$((180-i))/output/min_$((180-i)).out \
			  -p tleap/ala_tri.prmtop \
			  -c samp_$((180-i+3))/output/samp_$((180-i+3)).rst7 \
			  -r min_$((180-i))/output/min$((180-i)).rst7

		$PMEMD -O -i equil_$((180-i))/equil_$((180-i)).in -o equil_$((180-i))/output/equil_$((180-i)).out \
			  -p tleap/ala_tri.prmtop -c min_$((180-i))/output/min$((180-i)).rst7 \
			  -r equil_$((180-i))/output/equil_$((180-i)).rst7
			  
		$PMEMD -O -i samp_$((180-i))/samp_$((180-i)).in -o samp_$((180-i))/output/samp_$((180-i)).out \
			  -p tleap/ala_tri.prmtop -c equil_$((180-i))/output/equil_$((180-i)).rst7 \
			  -r samp_$((180-i))/output/samp_$((180-i)).rst7 \
			  -x samp_$((180-i))/output/samp_$((180-i)).nc
	fi
done
