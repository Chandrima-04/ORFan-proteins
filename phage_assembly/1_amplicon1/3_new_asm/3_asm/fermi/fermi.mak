FERMI=fermi
UNITIG_K=40
OVERLAP_K=48

all:fmdef.p2.mag.gz

# Construct the FM-index for raw sequences
fmdef.raw.fmd:../../2_qf/4_phix_rm/d9539_amplicon1_qf_pe.fq.gz ../../2_qf/4_phix_rm/d9539_amplicon1_qf_se.fq.gz
	(gzip -dc ../../2_qf/4_phix_rm/d9539_amplicon1_qf_pe.fq.gz; gzip -dc ../../2_qf/4_phix_rm/d9539_amplicon1_qf_se.fq.gz) | $(FERMI) ropebwt -a bcr -v3 -btNf fmdef.raw.tmp - > $@ 2> $@.log

# Error correction
fmdef.ec.fq.gz:fmdef.raw.fmd
	(gzip -dc ../../2_qf/4_phix_rm/d9539_amplicon1_qf_pe.fq.gz; gzip -dc ../../2_qf/4_phix_rm/d9539_amplicon1_qf_se.fq.gz) | $(FERMI) correct -t 16  $< - 2> $@.log | gzip -1 > $@

# Construct the FM-index for corrected sequences
fmdef.ec.fmd:fmdef.ec.fq.gz
	$(FERMI) fltuniq $< 2> fmdef.fltuniq.log | $(FERMI) ropebwt -a bcr -v3 -btf fmdef.ec.tmp - > $@ 2> $@.log

# Generate unitigs
fmdef.p0.mag.gz:fmdef.ec.fmd
	$(FERMI) unitig -t 16 -l $(UNITIG_K) $< 2> $@.log | gzip -1 > $@

fmdef.p1.mag.gz:fmdef.p0.mag.gz
	$(FERMI) clean $< 2> $@.log | gzip -1 > $@
fmdef.p2.mag.gz:fmdef.p1.mag.gz
	$(FERMI) clean -CAOFo $(OVERLAP_K) $< 2> $@.log | gzip -1 > $@

