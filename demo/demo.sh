#!/bin/bash

# Genescoring
echo "~~~ Genescoring demo"
../pascalx -g False -c [1] -p 1 ensemble.txt EUR.simulated out/g_out genescoring -sh False -cr 0 -cp 4 gwasA.tsv.gz

echo "~~~ Pathway scoring demo"
../pascalx -g False -c [1] -p 1 -pw pw_test.tsv ensemble.txt EUR.simulated out/gp_out genescoring -sh False -cr 0 -cp 4 gwasA.tsv.gz

# X-scoring (coherence)
# No gpu, chr 1, cpu 1 
echo "~~~ X-scoring demo (coherence)"
../pascalx  -g False -c [1] -p 1 ensemble.txt EUR.simulated out/xc_out xscoring -sh1 False -cr1 0 -cp1 4 -cb1 3 -ca11 1 -ca21 2 gwasA.tsv.gz -sh2 False -cr2 0 -cp2 4 -cb2 3 -ca12 1 -ca22 2 gwasB.tsv.gz

# X-scoring (anti-coherence)
# No gpu, chr 2, cpu 1
echo "~~~ X-scoring demo (anti-coherence)"
../pascalx  -g False -c [1] -p 1 ensemble.txt EUR.simulated out/xa_out xscoring -sh1 False -t True -cr1 0 -cp1 4 -cb1 3 -ca11 1 -ca21 2 gwasA.tsv.gz -sh2 False -cr2 0 -cp2 4 -cb2 3 -ca12 1 -ca22 2 gwasB.tsv.gz

# Ratio-scoring (coherence)
# No gpu, chr 1, cpu 1 
echo "~~~ Ratio-scoring demo (coherence)"
../pascalx  -g False -c [1] -p 1 ensemble.txt EUR.simulated out/rc_out xscoring -sh1 False -r True -cr1 0 -cp1 4 -cb1 3 -ca11 1 -ca21 2 gwasA.tsv.gz -sh2 False -cr2 0 -cp2 4 -cb2 3 -ca12 1 -ca22 2 gwasB.tsv.gz

# Ratio-scoring (coherence, flipped)
# No gpu, chr 1, cpu 1 
echo "~~~ Ratio-scoring demo (coherence, flipped)"
../pascalx  -g False -c [1] -p 1 ensemble.txt EUR.simulated out/rcf_out xscoring -sh1 False -r True -f True -cr1 0 -cp1 4 -cb1 3 -ca11 1 -ca21 2 gwasA.tsv.gz -sh2 False -cr2 0 -cp2 4 -cb2 3 -ca12 1 -ca22 2 gwasB.tsv.gz
