  j=10

 while [ $j -le 20 ]; do
#   cp sig_output_t.00"$j" sig_output
#   cp converged_mu_t.00"$j" converged_mu

   sed s/"TEMP"/"0.""$j""d0"/ input_file  > in.tmp
   qsub p-wave.qsub

   while [ "0" -eq "0" ]; do
     alive=`qstat -a | grep p-wave`
     if [ "$alive" = "" ]
     then
        break
     else
        sleep 60
     fi
   done

   mv output_file output_file_t."$j"
   cp sig_output sig_output_t."$j"
   cp converged_mu converged_mu_t."$j"

   j=`expr $j + 5`

 done

