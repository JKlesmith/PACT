# PACT - Protein Analysis and Classifier Toolkit<br />

This software is released under GNU General Public License 3<br/>
Additional license options available at http://license.umn.edu<br/>

Contact: Justin R. Klesmith through the issue tracker.<br/>
# Citation:<br/>
https://github.com/JKlesmith/PACT/<br/>
https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty1042/5258100<br/>
https://doi.org/10.1093/bioinformatics/bty1042<br/>

Improved mutant function prediction via PACT: Protein Analysis and Classifier Toolkit <br/>
Justin R Klesmith  Benjamin J Hackel<br/>
Bioinformatics, bty1042<br/>

# Update: 2022.1 - revision 1
A key change in the enrichment calculation has been made. Previously during software development we decided that a mutation (while not in the library design) that passed the mutational count filter would be considered in the total count. It is theoretically possible that a poorly made DNA library leading to a lot of rejected mutations could skew the raw enrichments. However, this shouldn't be an issue as we typically use normalized fitness values. A sample from the CD19 epitope 3B10 dataset showed this to be the case (where the log2 enrichments showed minor differences in the fraction (yet regressed to 1 using a linear fit) but the normalized Z-score were exacty the same). Therefore, the new default action of the enrichment module is to not consider mutations not desiged in the library in the total count. If the old method is desired a line stating "consider_rejected: True" in the config file under the enrichment heading can be added.<br />

A second change is addition of two new modes for strict counting. Setting strict_count_threshold under the enrichment section to <i>true</i> or <i>both</i> will enforce strict counting for the reference and selected with the thresholds (no virtual counts of 1 if the corresponding reference or selected count passes its threshold). The new addition is setting <i>ref-only</i> or <i>sel-only</i> to enable strict counting for the reference and selected respectively and virtual counts for the corresponding. Setting a value of false will allow existing virtual count behavior for both.

# Usage:<br />
python pact.py -c ./path/to/config_file.ini<br/>
Example config files are stored in ./pact/tests/

# Main Protocols:<br />
back_to_consensus<br />
classification_analysis<br />
function_filter<br />
fitness<br />
pact_vs_pact<br />
pact_vs_analysis<br />
sequence_homology<br />
shannon_entropy<br />
structure_analysis<br />
tools<br />

# Experimental Protocols:<br />
epitope_mapping<br />
enzyme_solubility<br />

# Release History:<br/>
Version 2018.6 - Initial Release<br/>
