unset LD_LIBRARY_PATH
python annotation.py --genome fna/Macaca_mulatta_new_genomic.fna.masked --model_path ./ANNEVO_model/ANNEVO_Mammalia.pt --output ./output/Macaca_mulatta_new.gff --threads 16
python annotation.py --genome fna/Microcebus_murinus_new_genomic.fna.masked --model_path ./ANNEVO_model/ANNEVO_Mammalia.pt --output ./output/Microcebus_murinus_new.gff --threads 16
python annotation.py --genome fna/Callithrix_jacchus_new_genomic.fna.masked --model_path ./ANNEVO_model/ANNEVO_Mammalia.pt --output ./output/Callithrix_jacchus_new.gff --threads 16 --genome_size_threshold 20*1024*1024
for species in Callithrix_jacchus Macaca_mulatta Microcebus_murinus; do /home/lulab2025/anaconda3/envs/primate/bin/gffread "output/${species}_new.gff" -g "fna/${species}_new_genomic.fna.masked" -y "output/${species}_new.pep.fa"; done

python annotation.py --genome fna_tba/Eulemur_mongoz_genomic.fna.masked --model_path ./ANNEVO_model/ANNEVO_Mammalia.pt \
--output ./output1/Eulemur_mongoz_genomic.fna.masked.gff \
--threads 16 \
--genome_size_threshold 10485760 \
--batch_size 16 \
--num_workers 2 \
> ./logs/Eulemur_mongoz_genomic.fna.masked.log 2>&1

LOG=logs/SPN_windowcap_defaultbp_$(date +%F_%H%M%S).log && \
set -o pipefail && \
python -u annotation.py \
  --genome fna_tba/SPN_genomic.fna.masked \
  --model_path ./ANNEVO_model/ANNEVO_Mammalia.pt \
  --output ./output1/SPN_genomic.fna.masked.gff \
  --threads 16 \
  --genome_size_threshold 104857600 \
  --batch_size 32 \
  --max_windows_per_chunk 8192 \
  --num_workers 2 \
  2>&1 | tee "$LOG"