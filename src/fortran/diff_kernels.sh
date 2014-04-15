for i in `ls *.f90`;do
  vimdiff $i ~/Projects/CloverLeaf/Code/CloverLeaf_ref/$i && read "Press [Enter] to continue..."
done
