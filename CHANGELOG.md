# Changelog

01/05/2022 - Fixed bug in the likelihood method R scripts which occurred when processing a single query sample.

06/2024 Version: 1.2.1
- added the "--no-BAQ" parameter to samtools mpileup to prevent base quality recalibration. This option is known to reduce reference bias.
- added contamin.py script which estimates Y-chromosome contamination by looking at mismatches at ISOGG SNP sites. See tutorial for details.

13/09/2024 Version: 1.2.2
- fixed paths issue with conda installation

16/09/2024 Version: 1.2.3
- fixed broken samtools when installing with conda (specify samtools=1.20)
- fixed pathphynder help command
- added --version command
- added log

07/06/2025 Version: 1.2.4
 - updated haplogroup nomenclature in newick and VCF files to match the one used on ISOGG list of variants
(BigTree.Y_GRCh37.V2.250606.nwk and BigTree.Y_GRCh37.V2.250606.vcf.gz)

21/07/2026 Version: 1.3.0
 - rewrote determineHG() haplogroup-calling logic: tree-walk based determination with ancestor-contradiction rejection, evidence-weighted ambiguity resolution, and absolute-count qualification for long branches
 - added a topology-aware note surfacing a viable deeper placement when it isn't reachable by name-matching alone
 - added a "Rejected" section showing why a downstream candidate was not adopted, and ancestor/descendant context labels in Downstream conflicts
 - added -x/--extra_reports flag to skip writing diagnostic-only reports (all_paths_report.txt, hg_in_tree_status.txt, .all_counts, snp_info_along_path) by default
 - cleaned up console output (removed leftover debug prints, suppressed stray auto-prints)
 - restructured documentation into docs/installation.md, docs/workflow.md, and CHANGELOG.md; added scripts/update_workflow_options.sh to keep the CLI reference in sync with --help
