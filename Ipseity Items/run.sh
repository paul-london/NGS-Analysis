snakemake --cores 16 --keep-going --snakefile /home/ipseity-gnt-i9/ACMG/script/snakefile_strand_bias_final_plot --configfile ./config.yaml --rerun-incomplete
mv /home/ipseity-gnt-i9/ACMG/VEP_data/*_vep.vcf ./
python3 ../run_filter.py
#python3 ..depth_no_save.py
python3 ./depth_save.py
bash ../pseudo.sh
python3 ../APOE_geno.py > APOE_data.txt
python3 ../APOE_raw.py > APOE_raw_data.txt
python3 ../coverage_table.py
python3 ./depth_plot.py
python3 ./plot_pdf.py