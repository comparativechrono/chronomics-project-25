cd /home/claude/work
for i in $(seq 1 60); do
  python3 scripts/phaseC_variants.py 2>>out/phaseC.log && { echo "FINISHED pass $i"; break; }
  echo "pass $i incomplete, retrying"; sleep 2
done
echo "genes in csv: $(tail -n +2 out/promoter_snvs.csv | cut -d, -f1 | sort -u | wc -l)/642"
