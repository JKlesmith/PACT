# PACT - Protein Analysis and Classifier Toolkit<br />

This software is released under GNU General Public License 3<br/>
Additional license options available at http://license.umn.edu<br/>

Contact: Justin R. Klesmith at: jrk [at] umn.edu<br/>
# Citation:<br/>
https://github.com/JKlesmith/PACT/<br/>
https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty1042/5258100<br/>
https://doi.org/10.1093/bioinformatics/bty1042<br/>

Improved mutant function prediction via PACT: Protein Analysis and Classifier Toolkit <br/>
Justin R Klesmith  Benjamin J Hackel<br/>
Bioinformatics, bty1042<br/>

# Update: May 2019<br />
Moving most of the core to Numpy/Scipy.<br/>
Moving forward a new directory will be uploaded with scripts from associated papers that use PACT calculated data.<br/>
Suggest using PACT for dataset generation then use Spyder for individual analyses.<br/>

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
