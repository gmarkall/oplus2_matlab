for i in 64 128 256 512 1024; do
  echo "PARTITION SIZE = $i";
  for j in 32 64 128 256 512 1024; do
    echo "BLOCK SIZE = $j `./airfoil_cuda OP_BLOCK_SIZE=$j OP_PART_SIZE=$i | tail -7`";
    echo "++++++++++++++++++++++++++++++++++++++++++++++++"; 
  done
  echo "================================================"; 
  echo " "; 
done

