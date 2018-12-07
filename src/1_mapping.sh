./nextflow src/mapping.nf -c src/mapping.config -profile docker -resume -w ~/data/work/ --fasta "data/fasta/final_assembly.fasta" --fastq "data/fastq/*_R{1,2}.fastq.gz" --contig '["Backbone_1017", "Backbone_1061", "Backbone_1061", "Backbone_1092", "Backbone_1094", "Backbone_1097", "Backbone_1097", "Backbone_1112", "Backbone_1271", "Backbone_1277", "Backbone_1277", "Backbone_1277", "Backbone_1301", "Backbone_1318", "Backbone_1338", "Backbone_1482", "Backbone_1493", "Backbone_1579", "Backbone_162", "Backbone_1702", "Backbone_177", "Backbone_210", "Backbone_286", "Backbone_476", "Backbone_796", "Backbone_800", "Backbone_802", "Backbone_901", "Backbone_975", "Backbone_1504", "Backbone_1689", "Backbone_656"]'
