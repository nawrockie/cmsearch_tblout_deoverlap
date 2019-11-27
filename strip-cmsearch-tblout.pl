while($line = <>) { 
  if($line !~ m/^\#/){ 
    chomp $line;
#NC_013790.1          -         LSU_trypano_mito     RF02546   hmm      127      330  2472190  2472399      +     -    6 0.24  11.6   10.7      0.24 ?   Methanobrevibacter ruminantium M1 chromosome, complete genome
    @el_A = split(/\s+/, $line);
    if(scalar(@el_A) >= 18) { 
      printf("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s -\n",
             $el_A[0], $el_A[1], $el_A[2], $el_A[3], $el_A[4], $el_A[5], $el_A[6], $el_A[7], $el_A[8], $el_A[9], $el_A[10], $el_A[11], $el_A[12], $el_A[13], $el_A[14], $el_A[15], $el_A[16]);
    }
  }
}

      
