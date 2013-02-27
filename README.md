DegradomeAnalyseR:
=======

DegradomeAnalyseR is copyright 2012 Kevin Murray, and is licensed under the GPLv3 License.  

See http://www.gnu.org/licenses/gpl-3.0.txt for the full license text


USAGE:
======
These scripts are used in a bash pipeline, e.g.:

    - scrips/sanitise_paresnip_output.py <raw_paresnip_output> > <sanitised_file>
    - paresnip/paresnip_filter.R --args <file1> [<file2> <file3> ... <fileN>] <prefix>
    - paresnip/paresnip_summarise.R --args <conserved_peaks_table> <prefix>
    - rnaseq/paresnip_downstream_bias.R --args <summarised_peaks_table> <target_length_table> <RNAseq_filtered_file> <prefix>
    - rnaseq/filter_bam.R --args <bam_file> <prefix>


File pipeline:
==============

Overall pipeline:
-----------------
    - RNAseq + Paresnip happen independently
    - then, paresnip_downstream_bias.R calculates downstream bias

Paresnip results:
-----------------
    - scripts/sanitise_paresnip.py
    - paresnip/paresnip_filter.R
    - paresnip/paresnip_summarise.R
    - rnaseq/paresnip_correlate_FA_SRA.R

RNAseq:
-------
    - Align to reference genome
    - If required, "samtools view -bS" to convert sam to bam
    - filter_bam.R generates target and alignment databases
