for i in `qstat | tail -n +3`;
  do
  if [ "${i:0:2}" = "29" ];
  then
    qdel $i
  fi
done
