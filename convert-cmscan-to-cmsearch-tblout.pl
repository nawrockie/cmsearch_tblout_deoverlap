while($line = <>) { 
  if($line !~ m/^\#/){ 
    chomp $line;
##idx target name            accession query name           accession clan name mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc olp anyidx afrct1 afrct2 winidx wfrct1 wfrct2 description of target
##--- ---------------------- --------- -------------------- --------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- --- ------ ------ ------ ------ ------ ------ ---------------------
#1    LSU_rRNA_archaea       RF02540   NC_013790.1          -         -         hmm       10     2982   762882   765854      +     -    6 0.49  11.2 2124.0         0  !   ^       -      -      -      -      -      - Archaeal large subunit ribosomal RNA
    @el_A = split(/\s+/, $line);
    if(scalar(@el_A) >= 27) { 
      printf("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s -\n",
             $el_A[3], $el_A[4], $el_A[1], $el_A[2], $el_A[6], $el_A[7], $el_A[8], $el_A[9], $el_A[10], $el_A[11], $el_A[12], $el_A[13], $el_A[14], $el_A[15], $el_A[16], $el_A[17], $el_A[18]);
    }
    else { 
      printf("nels: %d\n", scalar(@el_A));
    }
  }
}

      
