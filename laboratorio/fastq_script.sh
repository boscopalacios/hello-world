#!/bin/bash
for laboratorio in *.fastq
   do echo $laboratorio
   wc -l $laboratorio
echo terminado 
done

